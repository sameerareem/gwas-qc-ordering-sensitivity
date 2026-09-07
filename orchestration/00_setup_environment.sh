#!/usr/bin/env bash
# =============================================================================
#  00_setup_environment.sh
#  Environment setup, dependency check, and directory initialisation
#  Run this FIRST before any other pipeline step.
# =============================================================================

set -euo pipefail

echo "=================================================================="
echo "  GWAS Synthetic Data Pipeline — Environment Setup"
echo "=================================================================="

# ── Directory structure ──
DIRS=(
  "./gwas_synthetic_output"
  "./gwas_synthetic_output/logs"
  "./gwas_synthetic_output/merged"
  "./gwas_synthetic_output/chrX"
  "./gwas_synthetic_output/plink_chroms"
  "./gwas_synthetic_output/hapgen2_raw"
  "./gwas_synthetic_output/phenotypes"
  "./gwas_synthetic_output/replicates"
  "./reference_data/1000G"
)
for d in "${DIRS[@]}"; do
  mkdir -p "$d"
  echo "[DIR] Created: $d"
done

echo ""
echo "── Checking required binaries ──"

check_binary() {
  local name=$1
  local cmd=${2:-$1}
  if command -v "$cmd" &>/dev/null; then
    local ver
    ver=$("$cmd" --version 2>&1 | head -1 || echo "version unknown")
    echo "[OK]  $name: $ver"
  else
    echo "[MISSING] $name — Please install $name and ensure it is in PATH"
    echo "          → See README.md for installation instructions"
  fi
}

check_binary "PLINK2"    "plink2"
check_binary "HAPGEN2"   "hapgen2"
check_binary "bcftools"  "bcftools"
check_binary "tabix"     "tabix"
check_binary "wget"      "wget"
check_binary "Python3"   "python3"
check_binary "Rscript"   "Rscript"

echo ""
echo "── Checking Python packages ──"
python3 - <<'EOF'
required = ["numpy", "pandas", "scipy", "yaml", "json", "pathlib"]
for pkg in required:
    try:
        __import__(pkg if pkg != "yaml" else "yaml")
        print(f"[OK]  {pkg}")
    except ImportError:
        print(f"[MISSING] {pkg}  → pip install {pkg}")
EOF

echo ""
echo "── Checking R packages ──"
Rscript - <<'EOF'
required <- c("PhenotypeSimulator", "data.table", "optparse", "jsonlite")
for (pkg in required) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    ver <- packageVersion(pkg)
    cat(sprintf("[OK]  %s v%s\n", pkg, ver))
  } else {
    cat(sprintf("[MISSING] %s  → install.packages('%s')\n", pkg, pkg))
  }
}
EOF

echo ""
echo "── HAPGEN2 Installation Note ──"
echo "  HAPGEN2 is not on CRAN/PyPI. Download from:"
echo "  https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html"
echo "  Then: chmod +x hapgen2 && sudo mv hapgen2 /usr/local/bin/"
echo ""
echo "── 1000G Reference Files Note ──"
echo "  IMPUTE2-format haplotype files for HAPGEN2 are available at:"
echo "  https://mathgen.stats.ox.ac.uk/impute/data_download/1000GP_Phase3/"
echo "  These are ~20GB total. The pipeline downloads them automatically in Step 1."
echo ""
echo "── Quick Install Commands ──"
echo "  # Python packages"
echo "  pip install numpy pandas scipy pyyaml"
echo ""
echo "  # PLINK2"
echo "  wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64.zip"
echo "  unzip plink2_linux_x86_64.zip && sudo mv plink2 /usr/local/bin/"
echo ""
echo "  # bcftools + tabix (via conda)"
echo "  conda install -c bioconda bcftools tabix"
echo ""
echo "  # R packages"
echo '  Rscript -e "install.packages(c('"'"'PhenotypeSimulator'"'"','"'"'data.table'"'"','"'"'optparse'"'"','"'"'jsonlite'"'"'))"'
echo ""
echo "=================================================================="
echo "  Setup check complete. Review any [MISSING] items above."
echo "  Then run: python gwas_pipeline_master.py --config config.yaml"
echo "=================================================================="
