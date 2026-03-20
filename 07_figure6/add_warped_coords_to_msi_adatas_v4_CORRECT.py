#!/usr/bin/env python
"""
Add Warped Coordinates to MSI AnnData Files - Version 4 (CORRECT)

ROOT CAUSE IDENTIFIED:
- The MSI->DAPI registrar used the WARPED DAPI (already in H&E space, 11500x6000)
- This was copied from: out_feb2025_xenium/.../registered/morphology_focus_valis.ome.tiff
- So MSI coords warped through this registrar are ALREADY in H&E space!
- We do NOT need to chain through DAPI->H&E registration.

The correct coordinate chain is:
  MSI native pixels
  -> (pad 50px, transpose, flip)
  -> VALIS MSI input
  -> (MSI->DAPI registrar, but DAPI is actually warped-DAPI = H&E space)
  -> H&E space (11500 x 6000 cropped region)

For Xenium:
  Xenium microns
  -> DAPI pixels (divide by pixel size)
  -> (DAPI->H&E registrar)
  -> H&E space

Now both MSI and Xenium should be in the same H&E coordinate space!

Author: Generated for Chiba MDACC Project
Date: 2026-01-30
"""

import os
import numpy as np
import anndata as ad
from pathlib import Path
from valis import registration
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================

BASE_DIR = Path(r'T:\Analysis\116_MDACC_Chiba')
CHIBA_DIR = BASE_DIR / 'chiba'
NEW_CHIBA_DIR = BASE_DIR / 'New_Chiba_Jan26' / 'chiba-project-vm-jan26'

# VALIS registrar paths
MSI_TO_HE_REGISTRAR = (
    CHIBA_DIR / 'out_feb2025_msi' / 'script02b_registration_out_msi_to_dapi' /
    'results' / 'script02a_converted_images_for_valis_reg' / 'data' /
    'script02a_converted_images_for_valis_reg_registrar.pickle'
)  # Note: "DAPI" in this registrar is actually H&E-space DAPI!

DAPI_TO_HE_REGISTRAR = (
    CHIBA_DIR / 'out_feb2025_xenium' / 'script01b_registration_out_xeDapi_he__high_res_he_' /
    'results' / 'reg01_dapi_he' / 'data' / 'reg01_dapi_he_registrar.pickle'
)

# Xenium pixel size (microns per pixel in the morphology images)
XENIUM_PIXEL_SIZE_UM = 0.2125

# MSI slide names in registrar
MSI_SLIDE_NAMES = {
    'glycans': 'maldi_gycans_rep',
    'mets': 'mets_rep',
    'peptides': 'peptides_rep',
}

# Input/Output directories
MSI_INPUT_DIR = NEW_CHIBA_DIR / 'out_feb2025_msi'
XENIUM_DIR = CHIBA_DIR / 'adatas_with_rctd'
OUTPUT_DIR = NEW_CHIBA_DIR / 'out_feb2025_msi' / 'script03e_adatas_msi_aligned_CORRECT'

OVERWRITE = True

# ============================================================================
# Functions
# ============================================================================

def load_registrars():
    """Load VALIS registrars."""
    print("Loading VALIS registrars...")
    msi_to_he = registration.load_registrar(str(MSI_TO_HE_REGISTRAR))
    dapi_to_he = registration.load_registrar(str(DAPI_TO_HE_REGISTRAR))

    # Verify our understanding
    dapi_msi_reg = msi_to_he.get_slide('morphology_focus_valis')
    dapi_he_reg = dapi_to_he.get_slide('morphology_focus_valis')
    he_slide = dapi_to_he.get_slide('he_rotated_3dim')

    print(f"\nCoordinate Space Verification:")
    print(f"  MSI Registrar 'DAPI' (actually warped-DAPI): {dapi_msi_reg.slide_shape_rc}")
    print(f"  DAPI->H&E Registrar original DAPI: {dapi_he_reg.slide_shape_rc}")
    print(f"  DAPI->H&E Registrar H&E: {he_slide.slide_shape_rc}")

    # The MSI registrar's "DAPI" should match the H&E dimensions
    if np.allclose(dapi_msi_reg.slide_shape_rc, he_slide.slide_shape_rc):
        print(f"  [CONFIRMED] MSI 'DAPI' shape matches H&E shape -> MSI outputs are in H&E space!")
    else:
        print(f"  [WARNING] Shapes don't match exactly")

    return msi_to_he, dapi_to_he


def warp_msi_to_he(coords, modality, msi_to_he_reg):
    """
    Warp MSI coordinates to H&E space.

    The MSI registrar's 'DAPI' target is actually the DAPI warped to H&E space,
    so the output of this function is directly in H&E coordinates!
    """
    msi_slide = msi_to_he_reg.get_slide(MSI_SLIDE_NAMES[modality])
    # The target is called 'morphology_focus_valis' but it's actually in H&E space
    he_space_slide = msi_to_he_reg.get_slide('morphology_focus_valis')

    coords = np.asarray(coords, dtype=np.float64)
    coords_he = msi_slide.warp_xy_from_to(coords, to_slide_obj=he_space_slide)

    return coords_he


def warp_xenium_to_he(x_um, y_um, dapi_to_he_reg, pixel_size=XENIUM_PIXEL_SIZE_UM):
    """
    Warp Xenium coordinates (in microns) to H&E space.

    Steps:
    1. Convert Xenium microns to DAPI pixels
    2. Warp through DAPI->H&E registrar
    """
    # Convert microns to DAPI pixels
    x_dapi = x_um / pixel_size
    y_dapi = y_um / pixel_size

    coords_dapi = np.column_stack([x_dapi, y_dapi])

    # Warp from DAPI to H&E
    dapi_slide = dapi_to_he_reg.get_slide('morphology_focus_valis')
    he_slide = dapi_to_he_reg.get_slide('he_rotated_3dim')

    coords_he = dapi_slide.warp_xy_from_to(coords_dapi, to_slide_obj=he_slide)

    return coords_he


def process_sample(sample_name, modality, msi_to_he_reg, dapi_to_he_reg):
    """Process a single sample: warp both MSI and Xenium to H&E space for comparison."""

    # Load MSI data
    msi_subdir = f'script03b_adatas_msi_{modality}'
    msi_path = MSI_INPUT_DIR / msi_subdir / f'adata_{sample_name}.h5ad'

    if not msi_path.exists():
        print(f"  MSI file not found: {msi_path.name}")
        return None

    adata_msi = ad.read_h5ad(msi_path)

    # Load corresponding Xenium data
    xenium_path = XENIUM_DIR / f'adata_{sample_name}_with_rctd.h5ad'
    if not xenium_path.exists():
        print(f"  Xenium file not found: {xenium_path.name}")
        return None

    adata_xe = ad.read_h5ad(xenium_path)

    # Get MSI coordinates (original)
    if 'spatial' in adata_msi.obsm:
        msi_coords = adata_msi.obsm['spatial'].copy()
    else:
        msi_coords = np.column_stack([adata_msi.obs['x'].values, adata_msi.obs['y'].values])

    # Warp MSI coords to H&E space (directly, since target is warped-DAPI = H&E space)
    msi_coords_he = warp_msi_to_he(msi_coords, modality, msi_to_he_reg)

    # Get Xenium coordinates and warp to H&E space
    xe_x_um = adata_xe.obs['x_centroid'].values
    xe_y_um = adata_xe.obs['y_centroid'].values
    xe_coords_he = warp_xenium_to_he(xe_x_um, xe_y_um, dapi_to_he_reg)

    # Store MSI coordinates in H&E space
    adata_msi.obsm['spatial_he'] = msi_coords_he
    adata_msi.obs['x_he'] = msi_coords_he[:, 0]
    adata_msi.obs['y_he'] = msi_coords_he[:, 1]

    # Store Xenium coordinates in H&E space for reference
    adata_msi.uns['xenium_coords_he'] = {
        'x': xe_coords_he[:, 0],
        'y': xe_coords_he[:, 1],
        'n_cells': len(xe_coords_he),
        'x_range': [float(xe_coords_he[:, 0].min()), float(xe_coords_he[:, 0].max())],
        'y_range': [float(xe_coords_he[:, 1].min()), float(xe_coords_he[:, 1].max())],
    }

    # Store original Xenium coordinates for reference
    adata_msi.uns['xenium_coords_microns'] = {
        'x_range': [float(xe_x_um.min()), float(xe_x_um.max())],
        'y_range': [float(xe_y_um.min()), float(xe_y_um.max())],
    }

    adata_msi.uns['transform_info'] = {
        'msi_to_he': 'MSI -> warped-DAPI (=H&E space) via VALIS',
        'xenium_to_he': 'Xenium microns -> DAPI pixels -> H&E via VALIS',
        'xenium_pixel_size_um': XENIUM_PIXEL_SIZE_UM,
        'note': 'Both MSI and Xenium coords are now in H&E space (11500x6000)',
    }

    # Print comparison - these should now overlap!
    print(f"  MSI in H&E:    X=[{msi_coords_he[:,0].min():.1f}, {msi_coords_he[:,0].max():.1f}], Y=[{msi_coords_he[:,1].min():.1f}, {msi_coords_he[:,1].max():.1f}]")
    print(f"  Xenium in H&E: X=[{xe_coords_he[:,0].min():.1f}, {xe_coords_he[:,0].max():.1f}], Y=[{xe_coords_he[:,1].min():.1f}, {xe_coords_he[:,1].max():.1f}]")

    # Check for overlap
    msi_bbox = [msi_coords_he[:,0].min(), msi_coords_he[:,0].max(),
                msi_coords_he[:,1].min(), msi_coords_he[:,1].max()]
    xe_bbox = [xe_coords_he[:,0].min(), xe_coords_he[:,0].max(),
               xe_coords_he[:,1].min(), xe_coords_he[:,1].max()]

    x_overlap = max(0, min(msi_bbox[1], xe_bbox[1]) - max(msi_bbox[0], xe_bbox[0]))
    y_overlap = max(0, min(msi_bbox[3], xe_bbox[3]) - max(msi_bbox[2], xe_bbox[2]))

    if x_overlap > 0 and y_overlap > 0:
        print(f"  [OK] Bounding boxes OVERLAP! X_overlap={x_overlap:.1f}, Y_overlap={y_overlap:.1f}")
    else:
        print(f"  [WARNING] No overlap detected!")

    return adata_msi


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 70)
    print("MSI-Xenium Alignment - Version 4 (CORRECT)")
    print("=" * 70)
    print("""
Key insight: The MSI registrar's "DAPI" target is actually the DAPI image
that was already warped to H&E space (11500x6000). So MSI warped coords
are directly in H&E space - no need to chain through DAPI->H&E!
""")

    # Load registrars
    msi_to_he_reg, dapi_to_he_reg = load_registrars()

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput: {OUTPUT_DIR}")

    # Get all samples from input directory
    all_samples = set()
    for modality in ['glycans', 'mets', 'peptides']:
        msi_subdir = MSI_INPUT_DIR / f'script03b_adatas_msi_{modality}'
        if msi_subdir.exists():
            for f in msi_subdir.glob('adata_*.h5ad'):
                sample = f.stem.replace('adata_', '')
                all_samples.add(sample)

    all_samples = sorted(all_samples)
    print(f"\nFound {len(all_samples)} samples: {all_samples[:5]}...")

    for modality in ['glycans', 'mets', 'peptides']:
        print(f"\n{'='*70}")
        print(f"Processing {modality.upper()}")
        print("=" * 70)

        for sample in all_samples:
            print(f"\n  {sample}:")

            result = process_sample(sample, modality, msi_to_he_reg, dapi_to_he_reg)

            if result is not None:
                # Save
                out_path = OUTPUT_DIR / modality / f'adata_{sample}.h5ad'
                out_path.parent.mkdir(parents=True, exist_ok=True)
                result.write_h5ad(out_path)
                print(f"  Saved: {out_path.name}")

    print("\n" + "=" * 70)
    print("DONE - Check if coordinates now align correctly!")
    print("=" * 70)


if __name__ == '__main__':
    main()
