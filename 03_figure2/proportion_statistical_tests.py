#!/usr/bin/env python3
"""
Cell Type Proportion Statistical Tests
=======================================
Manuscript: "Integrated single cell spatial multi-omics landscape of WHO grades 2-4 gliomas"
Figures: 2h-j (core vs. edge within each tumor type) and Supplementary Figure 3 (across tumor types)

This script performs:
1. Paired Wilcoxon signed-rank tests comparing cell type proportions between
   core and edge regions within each tumor type (Figure 2h-j).
   Core and edge samples are paired by patient.

2. Kruskal-Wallis tests followed by pairwise Wilcoxon rank-sum tests comparing
   cell type proportions across tumor types (Supplementary Figure 3c-d).

All p-values are corrected for multiple comparisons using the Benjamini-Hochberg
(FDR) procedure within each comparison family.

Requirements:
    - pandas, numpy, scipy, statsmodels
    - RCTD output files (cell type proportions per sample)

Usage:
    python proportion_statistical_tests.py --rctd_dir <path_to_rctd_outputs> --output_dir <path_to_output>

Author: Trevor McKee / Astraea Bio
Date: 2026-03
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


# ──────────────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────────────

# Tumor type groupings — maps sample identifiers to tumor type and region
# Update these mappings to match your actual sample naming convention
TUMOR_TYPES = {
    'Oligo': 'IDHmut Oligodendroglioma',
    'AA': 'IDHmut Anaplastic Astrocytoma',
    'GBM': 'IDHwt Glioblastoma',
}

# Significance threshold after BH correction
ALPHA = 0.05


# ──────────────────────────────────────────────────────────────────────────────
# Data Loading
# ──────────────────────────────────────────────────────────────────────────────

def load_proportions(rctd_dir: str) -> pd.DataFrame:
    """
    Load RCTD cell type proportion data from the analysis output.

    Expected input: A CSV or set of CSVs with columns:
        - sample_id: Unique identifier for each TMA core
        - patient_id: Patient identifier (for pairing core/edge)
        - tumor_type: One of 'Oligo', 'AA', 'GBM'
        - region: 'Core' or 'Edge'
        - <cell_type_1>, <cell_type_2>, ...: Proportion columns (0-1 scale)

    If your data is in a different format, modify this function accordingly.

    Parameters
    ----------
    rctd_dir : str
        Path to directory containing RCTD proportion output files.

    Returns
    -------
    pd.DataFrame
        Tidy dataframe with proportions per sample.
    """
    rctd_path = Path(rctd_dir)

    # Try loading a single combined file first
    combined_file = rctd_path / 'cell_type_proportions_all_samples.csv'
    if combined_file.exists():
        df = pd.read_csv(combined_file)
        print(f"  Loaded {len(df)} samples from {combined_file}")
        return df

    # Otherwise, try loading individual per-sample files and combining
    csv_files = sorted(rctd_path.glob('*.csv'))
    if not csv_files:
        raise FileNotFoundError(
            f"No CSV files found in {rctd_dir}. "
            "Please provide RCTD proportion outputs."
        )

    frames = []
    for f in csv_files:
        try:
            frame = pd.read_csv(f)
            frames.append(frame)
        except Exception as e:
            warnings.warn(f"Could not read {f}: {e}")

    if not frames:
        raise ValueError("No valid proportion files found.")

    df = pd.concat(frames, ignore_index=True)
    print(f"  Loaded {len(df)} samples from {len(csv_files)} files")
    return df


def validate_data(df: pd.DataFrame) -> tuple:
    """
    Validate the proportion dataframe and extract cell type columns.

    Returns
    -------
    tuple of (pd.DataFrame, list)
        Validated dataframe and list of cell type column names.
    """
    required_cols = {'sample_id', 'patient_id', 'tumor_type', 'region'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    # Cell type columns are everything except the metadata columns
    cell_type_cols = [c for c in df.columns if c not in required_cols]
    if not cell_type_cols:
        raise ValueError("No cell type proportion columns found.")

    print(f"  Found {len(cell_type_cols)} cell types: {cell_type_cols}")
    print(f"  Tumor types: {df['tumor_type'].unique()}")
    print(f"  Regions: {df['region'].unique()}")
    print(f"  Patients: {df['patient_id'].nunique()}")

    return df, cell_type_cols


# ──────────────────────────────────────────────────────────────────────────────
# Test 1: Paired Core vs. Edge within each tumor type (Figure 2h-j)
# ──────────────────────────────────────────────────────────────────────────────

def test_core_vs_edge_paired(df: pd.DataFrame, cell_type_cols: list) -> pd.DataFrame:
    """
    Paired Wilcoxon signed-rank test comparing core vs. edge proportions
    within each tumor type.

    Pairing: Core and edge samples from the same patient are paired.

    Multiple testing correction: BH/FDR within each tumor type across
    all cell types tested.

    Parameters
    ----------
    df : pd.DataFrame
        Proportion data with columns: sample_id, patient_id, tumor_type,
        region, and cell type proportion columns.
    cell_type_cols : list
        Names of columns containing cell type proportions.

    Returns
    -------
    pd.DataFrame
        Results table with columns: tumor_type, cell_type, n_pairs,
        median_core, median_edge, median_diff, statistic, p_value,
        p_value_adj_BH, significant_BH, direction.
    """
    results = []

    for tumor_type in sorted(df['tumor_type'].unique()):
        tumor_df = df[df['tumor_type'] == tumor_type].copy()

        # Get patients with BOTH core and edge
        core_patients = set(
            tumor_df[tumor_df['region'] == 'Core']['patient_id']
        )
        edge_patients = set(
            tumor_df[tumor_df['region'] == 'Edge']['patient_id']
        )
        paired_patients = sorted(core_patients & edge_patients)

        if len(paired_patients) < 3:
            warnings.warn(
                f"  {tumor_type}: Only {len(paired_patients)} paired patients "
                f"(need >= 3 for Wilcoxon signed-rank). Skipping."
            )
            continue

        print(f"\n  {tumor_type}: {len(paired_patients)} paired patients")

        for ct in cell_type_cols:
            core_vals = []
            edge_vals = []

            for pid in paired_patients:
                c = tumor_df[
                    (tumor_df['patient_id'] == pid) &
                    (tumor_df['region'] == 'Core')
                ][ct]
                e = tumor_df[
                    (tumor_df['patient_id'] == pid) &
                    (tumor_df['region'] == 'Edge')
                ][ct]

                if len(c) == 1 and len(e) == 1:
                    core_vals.append(c.values[0])
                    edge_vals.append(e.values[0])

            core_arr = np.array(core_vals)
            edge_arr = np.array(edge_vals)
            diff = edge_arr - core_arr

            n_pairs = len(core_arr)

            if n_pairs < 3:
                continue

            # Wilcoxon signed-rank test (two-sided)
            # Handle case where all differences are zero
            if np.all(diff == 0):
                stat, pval = np.nan, 1.0
            else:
                try:
                    stat, pval = stats.wilcoxon(
                        core_arr, edge_arr, alternative='two-sided'
                    )
                except ValueError:
                    # All differences are zero or sample too small
                    stat, pval = np.nan, 1.0

            direction = 'Edge > Core' if np.median(diff) > 0 else 'Core > Edge'
            if np.median(diff) == 0:
                direction = 'No difference'

            results.append({
                'tumor_type': tumor_type,
                'cell_type': ct,
                'n_pairs': n_pairs,
                'median_core': np.median(core_arr),
                'median_edge': np.median(edge_arr),
                'median_diff': np.median(diff),
                'statistic': stat,
                'p_value': pval,
                'direction': direction,
            })

    if not results:
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # BH correction within each tumor type
    corrected_pvals = []
    for tumor_type in results_df['tumor_type'].unique():
        mask = results_df['tumor_type'] == tumor_type
        pvals = results_df.loc[mask, 'p_value'].values
        _, padj, _, _ = multipletests(pvals, alpha=ALPHA, method='fdr_bh')
        corrected_pvals.extend(padj)

    results_df['p_value_adj_BH'] = corrected_pvals
    results_df['significant_BH'] = results_df['p_value_adj_BH'] < ALPHA

    return results_df


# ──────────────────────────────────────────────────────────────────────────────
# Test 2: Across tumor types (Supplementary Figure 3)
# ──────────────────────────────────────────────────────────────────────────────

def test_across_tumor_types(df: pd.DataFrame, cell_type_cols: list) -> tuple:
    """
    Kruskal-Wallis test followed by pairwise Wilcoxon rank-sum tests
    comparing cell type proportions across tumor types, separately
    for core and edge regions.

    Multiple testing correction: BH/FDR across all cell types within
    each region (core or edge).

    Parameters
    ----------
    df : pd.DataFrame
        Proportion data.
    cell_type_cols : list
        Names of cell type proportion columns.

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame)
        (kruskal_results, pairwise_results)
    """
    kruskal_results = []
    pairwise_results = []

    tumor_types = sorted(df['tumor_type'].unique())

    for region in ['Core', 'Edge']:
        region_df = df[df['region'] == region].copy()

        print(f"\n  Region: {region}")
        for ct in cell_type_cols:
            groups = []
            group_labels = []
            for tt in tumor_types:
                vals = region_df[region_df['tumor_type'] == tt][ct].dropna().values
                if len(vals) > 0:
                    groups.append(vals)
                    group_labels.append(tt)

            if len(groups) < 2:
                continue

            # Kruskal-Wallis
            if len(groups) >= 2 and all(len(g) >= 1 for g in groups):
                try:
                    kw_stat, kw_pval = stats.kruskal(*groups)
                except ValueError:
                    kw_stat, kw_pval = np.nan, 1.0
            else:
                kw_stat, kw_pval = np.nan, 1.0

            kruskal_results.append({
                'region': region,
                'cell_type': ct,
                'n_groups': len(groups),
                'group_sizes': str([len(g) for g in groups]),
                'kw_statistic': kw_stat,
                'kw_p_value': kw_pval,
            })

            # Pairwise Wilcoxon rank-sum tests
            for i in range(len(groups)):
                for j in range(i + 1, len(groups)):
                    try:
                        rs_stat, rs_pval = stats.mannwhitneyu(
                            groups[i], groups[j], alternative='two-sided'
                        )
                    except ValueError:
                        rs_stat, rs_pval = np.nan, 1.0

                    pairwise_results.append({
                        'region': region,
                        'cell_type': ct,
                        'group_1': group_labels[i],
                        'group_2': group_labels[j],
                        'n_1': len(groups[i]),
                        'n_2': len(groups[j]),
                        'median_1': np.median(groups[i]),
                        'median_2': np.median(groups[j]),
                        'U_statistic': rs_stat,
                        'p_value': rs_pval,
                    })

    kw_df = pd.DataFrame(kruskal_results)
    pw_df = pd.DataFrame(pairwise_results)

    # BH correction within each region for Kruskal-Wallis
    if len(kw_df) > 0:
        for region in kw_df['region'].unique():
            mask = kw_df['region'] == region
            pvals = kw_df.loc[mask, 'kw_p_value'].values
            _, padj, _, _ = multipletests(pvals, alpha=ALPHA, method='fdr_bh')
            kw_df.loc[mask, 'kw_p_value_adj_BH'] = padj
        kw_df['kw_significant_BH'] = kw_df['kw_p_value_adj_BH'] < ALPHA

    # BH correction within each region for pairwise tests
    if len(pw_df) > 0:
        for region in pw_df['region'].unique():
            mask = pw_df['region'] == region
            pvals = pw_df.loc[mask, 'p_value'].values
            _, padj, _, _ = multipletests(pvals, alpha=ALPHA, method='fdr_bh')
            pw_df.loc[mask, 'p_value_adj_BH'] = padj
        pw_df['significant_BH'] = pw_df['p_value_adj_BH'] < ALPHA

    return kw_df, pw_df


# ──────────────────────────────────────────────────────────────────────────────
# Output
# ──────────────────────────────────────────────────────────────────────────────

def save_results(
    paired_df: pd.DataFrame,
    kruskal_df: pd.DataFrame,
    pairwise_df: pd.DataFrame,
    output_dir: str,
):
    """Save all results to CSV files."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    if len(paired_df) > 0:
        path = out / 'fig2hj_paired_core_vs_edge_proportions.csv'
        paired_df.to_csv(path, index=False, float_format='%.6g')
        print(f"\n  Saved: {path}")
        sig = paired_df[paired_df['significant_BH']]
        print(f"    {len(sig)} / {len(paired_df)} significant (BH adj. p < {ALPHA})")

    if len(kruskal_df) > 0:
        path = out / 'suppfig3_kruskal_wallis_across_tumors.csv'
        kruskal_df.to_csv(path, index=False, float_format='%.6g')
        print(f"\n  Saved: {path}")

    if len(pairwise_df) > 0:
        path = out / 'suppfig3_pairwise_wilcoxon_across_tumors.csv'
        pairwise_df.to_csv(path, index=False, float_format='%.6g')
        print(f"\n  Saved: {path}")
        sig = pairwise_df[pairwise_df['significant_BH']]
        print(f"    {len(sig)} / {len(pairwise_df)} significant (BH adj. p < {ALPHA})")


def print_summary(paired_df: pd.DataFrame):
    """Print a human-readable summary of paired test results."""
    if len(paired_df) == 0:
        print("\n  No paired test results.")
        return

    print("\n" + "=" * 80)
    print("SUMMARY: Paired Core vs. Edge Proportion Comparisons (Figure 2h-j)")
    print("=" * 80)

    for tumor_type in paired_df['tumor_type'].unique():
        subset = paired_df[paired_df['tumor_type'] == tumor_type]
        sig = subset[subset['significant_BH']]
        print(f"\n  {tumor_type}:")
        print(f"    Total cell types tested: {len(subset)}")
        print(f"    Significant (BH p < {ALPHA}): {len(sig)}")

        if len(sig) > 0:
            for _, row in sig.iterrows():
                print(
                    f"      {row['cell_type']}: "
                    f"median diff = {row['median_diff']:.4f}, "
                    f"p_adj = {row['p_value_adj_BH']:.2e}, "
                    f"{row['direction']}"
                )


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Cell type proportion statistical tests for glioma multi-omics paper'
    )
    parser.add_argument(
        '--rctd_dir',
        type=str,
        required=True,
        help='Directory containing RCTD cell type proportion CSV files',
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        default='./proportion_test_results',
        help='Directory for output CSV files (default: ./proportion_test_results)',
    )
    args = parser.parse_args()

    print("=" * 80)
    print("Cell Type Proportion Statistical Tests")
    print("=" * 80)

    # Load data
    print("\n[1/4] Loading proportion data...")
    df = load_proportions(args.rctd_dir)
    df, cell_type_cols = validate_data(df)

    # Test 1: Paired core vs. edge
    print("\n[2/4] Paired Wilcoxon signed-rank tests (Core vs. Edge)...")
    paired_df = test_core_vs_edge_paired(df, cell_type_cols)

    # Test 2: Across tumor types
    print("\n[3/4] Kruskal-Wallis + pairwise Wilcoxon rank-sum (across tumor types)...")
    kruskal_df, pairwise_df = test_across_tumor_types(df, cell_type_cols)

    # Save
    print("\n[4/4] Saving results...")
    save_results(paired_df, kruskal_df, pairwise_df, args.output_dir)

    # Summary
    print_summary(paired_df)

    print("\n" + "=" * 80)
    print("Done.")
    print("=" * 80)


if __name__ == '__main__':
    main()
