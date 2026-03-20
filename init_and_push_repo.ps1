# ============================================================================
# Initialize and push the glioma multi-omics code repository (PowerShell)
#
# Run from the repository root:
#   cd "T:\Analysis\116_MDACC_Chiba\Integrated_Spatial_Multiomics_Glioma_Repository"
#   .\init_and_push_repo.ps1
#
# Prerequisites:
#   1. Git installed and on PATH
#   2. GitHub credentials configured (gh auth login, or SSH key, or credential manager)
#   3. A new empty repo created on GitHub:
#      https://github.com/AstraeaBio/chiba-glioma-multiomics
# ============================================================================

# ---- CONFIGURATION - UPDATE THIS ----
$RemoteURL = "https://github.com/AstraeaBio/chiba_glioma_multiomics.git"
# For SSH: $RemoteURL = "git@github.com:AstraeaBio/chiba-glioma-multiomics.git"

# ---- Check git is available ----
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    Write-Error "Git is not installed or not on PATH. Install from https://git-scm.com/"
    exit 1
}

# ---- Initialize ----
Write-Host "Initializing git repository..." -ForegroundColor Cyan
git init
git branch -M main

# ---- Create .gitignore ----
$gitignore = @"
# Python
__pycache__/
*.py[cod]
*`$py.class
*.so
*.egg-info/
dist/
build/
.eggs/

# Jupyter checkpoints
.ipynb_checkpoints/

# R
.Rhistory
.Rdata
.RData
.Rproj.user/

# Environment
.env
*.env

# OS files
.DS_Store
Thumbs.db
desktop.ini

# Data files (too large for git)
*.h5ad
*.h5
*.ome.tiff
*.tiff
*.tif
*.imzML
*.ibd
*.csv.gz
*.pkl
*.pickle
*.rds
*.RData

# Output images (keep scripts, not generated figures)
*.png
*.jpg
*.jpeg
*.pdf
*.tiff
*.svg

# Large model files
*.pt
*.pth
*.ckpt

# IDE
.vscode/
.idea/
*.swp
*.swo

# Claude/AI assistant files
.claude/
CLAUDE.MD
"@

$gitignore | Out-File -FilePath ".gitignore" -Encoding utf8NoBOM
Write-Host ".gitignore created" -ForegroundColor Green

# ---- Create LICENSE ----
$license = @"
MIT License

Copyright (c) 2026 Astraea Bio LLC and The University of Texas MD Anderson Cancer Center

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"@

$license | Out-File -FilePath "LICENSE" -Encoding utf8NoBOM
Write-Host "LICENSE created" -ForegroundColor Green

# ---- Check for large files ----
Write-Host "`nChecking for large files (>10MB)..." -ForegroundColor Cyan
$largeFiles = Get-ChildItem -Recurse -File | Where-Object { $_.Length -gt 10MB -and $_.FullName -notlike "*\.git\*" }
if ($largeFiles) {
    Write-Host "  WARNING: The following files are >10MB and may need to be excluded:" -ForegroundColor Yellow
    foreach ($f in $largeFiles) {
        $sizeMB = [math]::Round($f.Length / 1MB, 1)
        Write-Host "    $($f.FullName) ($sizeMB MB)" -ForegroundColor Yellow
    }
    Write-Host ""
    $continue = Read-Host "Continue anyway? (y/n)"
    if ($continue -ne 'y') {
        Write-Host "Aborted. Add large files to .gitignore and try again."
        exit 0
    }
}

# ---- Stage everything ----
Write-Host "`nStaging files..." -ForegroundColor Cyan
git add -A

# ---- Show status ----
Write-Host "`nFiles staged for commit:" -ForegroundColor Cyan
git status --short
Write-Host ""

# ---- Commit ----
Write-Host "Creating initial commit..." -ForegroundColor Cyan
$commitMsg = @"
Initial commit: analysis code for Nature Medicine glioma multi-omics paper

Integrated single cell spatial multi-omics landscape of WHO grades 2-4 gliomas
identifies locoregional metabolomic surrogates of glioma transcriptional cellular states.

Ma Y, Ayyadhury S, Singh S, et al. Nature Medicine (2026).

Repository contains:
- Xenium spatial transcriptomics preprocessing and analysis (Python/scanpy)
- RCTD cell type deconvolution pipeline (R/spacexr)
- SCIMAP spatial cell-cell interaction analysis (Python)
- TGFb pathway analysis with statistical comparisons (Python)
- MSI-Xenium co-registration and metabolite integration (Python/VALIS)
- IMC segmentation and analysis (Python)
- Ki-67 metabolite enrichment analysis (Python)
- Cell type proportion statistical tests (Python)
- Environment configuration files

Note: Survival analysis scripts (UNT Denton) and MSPen validation
scripts (Baylor) to be added when received from collaborators.
"@

git commit -m $commitMsg

# ---- Add remote and push ----
Write-Host "`nAdding remote origin..." -ForegroundColor Cyan
git remote add origin $RemoteURL

Write-Host "`nPushing to GitHub..." -ForegroundColor Cyan
git push -u origin main

Write-Host ""
Write-Host "============================================" -ForegroundColor Green
Write-Host "Done! Repository pushed to:" -ForegroundColor Green
Write-Host "  $RemoteURL" -ForegroundColor White
Write-Host "============================================" -ForegroundColor Green
Write-Host ""
Write-Host "Next steps:" -ForegroundColor Cyan
Write-Host "  1. Verify the repo at the GitHub URL"
Write-Host "  2. Add a description and topics on GitHub"
Write-Host "  3. Update the manuscript Code Availability statement with the URL"
Write-Host "  4. When collaborator scripts arrive:"
Write-Host "       git add <file>"
Write-Host "       git commit -m 'Add survival analysis scripts from UNT'"
Write-Host "       git push"
Write-Host "  5. For the final submission, create a Zenodo DOI"
