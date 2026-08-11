#!/usr/bin/env python3
"""
================================================================================
 GWAS Synthetic Data Pipeline - Master Orchestration Script
 Study: Effect of QC Step Order on GWAS Phenotype Associations

 Pipeline Steps:
   1. Download 1000G reference haplotypes (autosomes + chrX)
   2. Simulate autosomal genotypes via msprime per chromosome
      (replaces HAPGEN2 which is incompatible with modern Linux kernels)
   3. Subsample & process real chrX data (1000G)
   4. Merge autosomes + chrX into unified PLINK dataset
   5. Simulate phenotypes (PhenotypeSimulator via R)
   6. Generate N datasets with randomised sex distortion (0-5%)
   7. Produce ground-truth causal-SNP and covariate files

 Usage:
   python gwas_pipeline_master.py --config config.yaml [--step STEP] [--dataset-id N]

 Requirements:
   - msprime >= 1.0, stdpopsim  (pip install msprime stdpopsim)
   - PLINK2  (binary in PATH or configured in config.yaml)
   - bcftools, tabix, bgzip (standard bioinformatics tools)
   - Python  >= 3.8  with: numpy, pandas, scipy, pyyaml
   - R       >= 4.0  with: PhenotypeSimulator, data.table

 Performance notes (all controlled via config.yaml):
   - Chromosomes simulated in parallel (n_parallel_chroms: 22)
   - Only 50MB subregion simulated per chromosome (region_length_bp)
   - Recombination map cropped to match subregion length (fixes ValueError)
   - Tree sequence thinned to target SNP count before VCF write
   - VCF streamed directly to bgzip — no uncompressed file written to disk
   - VCF deleted after PLINK conversion to save disk space
   - Output format: .bed/.bim/.fam (PLINK1 binary, standard for GWAS)

 Author: Generated for dissertation - QC Order Effects in GWAS
================================================================================
"""

import os
import sys
import math
import subprocess
import argparse
import logging
import yaml
import json
import random
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Optional
from multiprocessing import Pool
import traceback

# ─────────────────────────────────────────────────────────────────────────────
#  Setup
# ─────────────────────────────────────────────────────────────────────────────

def setup_logging(log_dir: str, dataset_id: Optional[int] = None) -> logging.Logger:
    """Initialise file + console logging."""
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    tag = f"_dataset{dataset_id}" if dataset_id is not None else ""
    log_file = os.path.join(log_dir, f"pipeline{tag}_{datetime.now():%Y%m%d_%H%M%S}.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger("gwas_pipeline")


def setup_worker_logging(log_dir: str, chrom: int) -> logging.Logger:
    """Per-chromosome logger for parallel worker processes."""
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    log_file = os.path.join(log_dir, f"chr{chrom}_{datetime.now():%Y%m%d_%H%M%S}.log")
    logger = logging.getLogger(f"chr{chrom}")
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        fh = logging.FileHandler(log_file)
        ch = logging.StreamHandler(sys.stdout)
        fmt = logging.Formatter(f"%(asctime)s [chr{chrom}] [%(levelname)s] %(message)s")
        fh.setFormatter(fmt)
        ch.setFormatter(fmt)
        logger.addHandler(fh)
        logger.addHandler(ch)
    return logger


def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def run_cmd(cmd: str, logger: logging.Logger, check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command, log it, and raise on failure."""
    logger.info(f"CMD: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.stdout and result.stdout.strip():
        logger.info(result.stdout.strip())   # INFO so PLINK messages are always visible
    if result.stderr and result.stderr.strip():
        logger.info(result.stderr.strip())   # INFO so warnings are never swallowed
    if check and result.returncode != 0:
        logger.error(f"Command failed (exit {result.returncode}):\n{result.stderr}")
        raise RuntimeError(f"Command failed: {cmd}")
    return result


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 1: Download 1000G Reference Data
# ─────────────────────────────────────────────────────────────────────────────

def download_1000g_reference(cfg: dict, logger: logging.Logger):
    """
    Download phased haplotype + legend files from 1000G Phase 3 (SHAPEIT2 format)
    and chrX VCF for real sex chromosome data.
    """
    ref_dir = Path(cfg["paths"]["reference_dir"])
    ref_dir.mkdir(parents=True, exist_ok=True)

    hapgen_ref_base = (
        "https://mathgen.stats.ox.ac.uk/impute/data_download/1000GP_Phase3"
    )

    logger.info("Downloading autosomal 1000G haplotype reference files ...")
    for chrom in cfg["reference"]["chromosomes"]:
        hap_file = ref_dir / f"1000GP_Phase3_chr{chrom}.hap.gz"
        leg_file = ref_dir / f"1000GP_Phase3_chr{chrom}.legend.gz"
        if not hap_file.exists():
            run_cmd(f"wget -q -P {ref_dir} {hapgen_ref_base}/1000GP_Phase3_chr{chrom}.hap.gz", logger)
        if not leg_file.exists():
            run_cmd(f"wget -q -P {ref_dir} {hapgen_ref_base}/1000GP_Phase3_chr{chrom}.legend.gz", logger)

    samp_file = ref_dir / "1000GP_Phase3.sample"
    if not samp_file.exists():
        run_cmd(f"wget -q -P {ref_dir} {hapgen_ref_base}/1000GP_Phase3.sample", logger)

    logger.info("Downloading chrX VCF from 1000G Phase 3 ...")
    chrX_vcf = ref_dir / "ALL.chrX.phase3.vcf.gz"
    chrX_tbi = ref_dir / "ALL.chrX.phase3.vcf.gz.tbi"
    tabix_bin = cfg["paths"].get("tabix_bin", "tabix")
    if not chrX_vcf.exists():
        url = (
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
            "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
        )
        # Use -O to save directly to target filename — avoids wildcard mv bug
        run_cmd(f"wget -q -O {chrX_vcf} {url}", logger)
    if not chrX_tbi.exists():
        run_cmd(f"{tabix_bin} -p vcf {chrX_vcf}", logger)

    logger.info("Reference download complete.")


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 2: Simulate Autosomal Genotypes via msprime
#  All four performance/space fixes applied:
#    1. Parallel chromosome simulation
#    2. 50MB region per chromosome instead of full chromosome
#    3. Thin tree sequence to target SNP count before writing VCF
#    4. Stream VCF directly to bgzip + delete after PLINK conversion
# ─────────────────────────────────────────────────────────────────────────────

def _simulate_one_chrom(args: tuple) -> str:
    """
    Worker function: simulate one chromosome with msprime and convert to PLINK2.
    Runs in a separate process when parallelised.
    Returns chrom_str on success, raises on failure.
    """
    chrom, cfg, work_dir_str, log_dir = args
    work_dir  = Path(work_dir_str)
    chrom_str = f"chr{chrom}"
    logger    = setup_worker_logging(log_dir, chrom)

    try:
        import msprime
        import stdpopsim
    except ImportError:
        raise RuntimeError("pip install msprime stdpopsim")

    # ── Paths ──
    vcf_dir   = work_dir / "msprime_vcf"
    plink_dir = work_dir / "plink_chroms"
    vcf_dir.mkdir(exist_ok=True)
    plink_dir.mkdir(exist_ok=True)

    vcf_out   = vcf_dir  / f"{chrom_str}.vcf.gz"
    plink_out = plink_dir / chrom_str

    # ── Skip if already converted ──
    if (plink_dir / f"{chrom_str}.bed").exists():
        logger.info(f"{chrom_str}: plink output exists, skipping.")
        return chrom_str

    # ── Read config values ──
    n_indiv        = cfg["n_individuals"]
    seed           = cfg["sex_distortion"]["seed"]
    mu             = cfg["msprime"].get("mutation_rate", 1.44e-8)
    plink2         = cfg["paths"]["plink2_bin"]
    maf_prefilter    = cfg["msprime"].get("plink_maf_prefilter", 0.001)
    region_length    = cfg["msprime"].get("region_length_bp", 50_000_000)
    thin_enabled     = cfg["msprime"].get("thin_to_target_snps", True)
    bgzip_bin        = cfg["paths"].get("bgzip_bin", "bgzip")
    tabix_bin        = cfg["paths"].get("tabix_bin", "tabix")
    target_snps    = cfg["reference"]["snps_per_chrom"][chrom]
    delete_vcf     = cfg["msprime"].get("delete_vcf_after_conversion", True)
    dem_model_name = cfg["msprime"].get("demographic_model", "OutOfAfrica_3G09")
    eur_pop_idx    = cfg["msprime"].get("eur_population_index", 1)

    # ── Load stdpopsim species + EUR demographic model ──
    logger.info(f"{chrom_str}: loading stdpopsim model '{dem_model_name}' ...")
    species = stdpopsim.get_species("HomSap")
    model   = species.get_demographic_model(dem_model_name)

    # ── Get contig with HapMap recombination map ──
    contig = species.get_contig(
        chrom_str,
        genetic_map="HapMapII_GRCh37",
    )

    # ── FIX 2: Use 50MB region instead of full chromosome ──
    # The recombination map must be cropped to exactly sim_length —
    # passing a shorter sequence_length with the full-length map raises:
    # ValueError: sequence_length must match the rate map length.
    sim_length = min(region_length, contig.length)

    full_map    = contig.recombination_map
    # Keep only breakpoints that fall within the subregion,
    # then cap the final position at sim_length exactly.
    crop_pos    = [p for p in full_map.left if p < sim_length] + [sim_length]
    crop_rates  = [full_map.rate[i]
                   for i, p in enumerate(full_map.left)
                   if p < sim_length]
    cropped_map = msprime.RateMap(position=crop_pos, rate=crop_rates)

    logger.info(
        f"{chrom_str}: simulating {n_indiv} individuals over "
        f"{sim_length/1e6:.0f}MB (full={contig.length/1e6:.0f}MB) ..."
    )

    # ── Run ancestry simulation ──
    ts = msprime.sim_ancestry(
        samples            = {eur_pop_idx: n_indiv},   # EUR population only
        demography         = model.model,
        recombination_rate = cropped_map,              # cropped map matches sim_length
        sequence_length    = sim_length,               # FIX 2: 50MB not full chr
        random_seed        = seed + chrom,
        ploidy             = 2,
        record_provenance  = False,                    # suppress 8MB provenance warning
    )

    # ── Overlay mutations ──
    # Note: num_threads was removed in newer msprime versions; threading
    # is handled internally and does not need to be specified.
    mts = msprime.sim_mutations(
        ts,
        rate        = mu,
        random_seed = seed + chrom + 1000,
        keep        = False,
    )

    logger.info(
        f"{chrom_str}: {mts.num_samples} haplotypes, "
        f"{mts.num_mutations} mutations before thinning."
    )

    # ── FIX 3: Thin tree sequence to target SNP count before writing VCF ──
    if thin_enabled and mts.num_sites > target_snps:
        all_sites  = list(mts.sites())
        step       = len(all_sites) // target_snps
        keep_idx   = set(range(0, len(all_sites), step)[:target_snps])
        drop_idx   = [i for i in range(mts.num_sites) if i not in keep_idx]
        mts        = mts.delete_sites(drop_idx)
        logger.info(
            f"{chrom_str}: thinned to {mts.num_sites} SNPs "
            f"(target={target_snps}, evenly spaced to preserve LD)."
        )

    # ── Assign SYNTH_ sample names before writing VCF ──
    # msprime default sample IDs are tsk_0, tsk_1, ... which are not unique
    # across chromosomes and do not match the SYNTH_ IDs used elsewhere.
    # We assign SYNTH_0000001 ... SYNTH_N names so all files share one ID system.
    n_samples    = mts.num_individuals
    synth_names  = [f"SYNTH_{i+1:07d}" for i in range(n_samples)]

    # ── Stream VCF directly to bgzip — never write uncompressed file ──
    # write_vcf() requires a text-mode stream but bgzip.stdin is binary.
    # Wrapping with TextIOWrapper bridges the str/bytes mismatch.
    import io
    logger.info(f"{chrom_str}: writing compressed VCF -> {vcf_out} ...")
    bgzip_proc = subprocess.Popen(
        [bgzip_bin, "-c"],              # use full conda path from config
        stdin  = subprocess.PIPE,
        stdout = open(vcf_out, "wb"),
        stderr = subprocess.PIPE,
    )
    text_stdin = io.TextIOWrapper(bgzip_proc.stdin, encoding="utf-8")
    mts.write_vcf(
        text_stdin,
        contig_id        = str(chrom),
        individual_names = synth_names,   # SYNTH_0000001 ... SYNTH_N
    )
    text_stdin.flush()
    text_stdin.detach()       # detach before close so stdin is not double-closed
    bgzip_proc.stdin.close()
    bgzip_proc.wait()
    if bgzip_proc.returncode != 0:
        raise RuntimeError(
            f"bgzip failed: {bgzip_proc.stderr.read().decode()}"
        )

    # Index the compressed VCF
    run_cmd(f"{tabix_bin} -p vcf {vcf_out}", logger)

    # ── Convert VCF → PLINK1 .bed/.bim/.fam ──
    # --set-all-var-ids forces unique chr:pos:ref:alt IDs for ALL variants.
    # --set-missing-var-ids is not sufficient because PLINK2 auto-assigns
    # sequential numeric IDs (6, 9, 10...) which are NOT unique across
    # chromosomes and cause PLINK1.9 merge to fail with
    # "Multiple chromosomes seen for variant" errors.
    cmd = (
        f"{plink2} "
        f"--vcf {vcf_out} "
        f"--chr {chrom} "
        f"--make-bed "              # outputs .bed/.bim/.fam
        f"--out {plink_out} "
        f"--maf {maf_prefilter} "   # loose pre-filter; real MAF done at QC stage
        f"--geno 0.05 "
        f"--max-alleles 2 "
        f"--snps-only "
        f"--set-all-var-ids @:#:$r:$a "            # force unique chr:pos:ref:alt IDs for ALL variants
        f"--new-id-max-allele-len 10 missing "      # truncate long alleles
    )
    run_cmd(cmd, logger)
    logger.info(f"{chrom_str}: PLINK2 conversion complete → {plink_out}")

    # ── Fix FID=0: set FID=IID in autosomal .fam ──
    # PLINK2 sets FID=0 for msprime VCFs (no family structure in VCF).
    # chrX blocks use FID=SYNTH_ so autosomes must match to avoid
    # PLINK1.9 treating (FID=0,IID=SYNTH_X) and (FID=SYNTH_X,IID=SYNTH_X)
    # as two different individuals during the final merge.
    fam_path = Path(str(plink_out) + ".fam")
    fam_lines = fam_path.read_text().splitlines()
    fixed = []
    needs_fix = False
    for line in fam_lines:
        parts = line.split()
        if parts[0] != parts[1]:  # FID != IID
            parts[0] = parts[1]   # set FID = IID
            needs_fix = True
        fixed.append("\t".join(parts))
    if needs_fix:
        fam_path.write_text("\n".join(fixed) + "\n")
        logger.info(f"{chrom_str}: FID set to IID in .fam")

    # ── FIX 4b: Delete VCF after conversion to free disk space ──
    if delete_vcf and vcf_out.exists():
        vcf_out.unlink()
        tbi = Path(str(vcf_out) + ".tbi")
        if tbi.exists():
            tbi.unlink()
        logger.info(f"{chrom_str}: VCF deleted to free disk space.")

    return chrom_str


def simulate_autosomes_hapgen2(cfg: dict, work_dir: Path, logger: logging.Logger):
    """
    Simulate autosomal genotypes for N individuals using msprime + stdpopsim.
    All 22 chromosomes run in parallel (FIX 1) using multiprocessing.Pool.
    Each worker applies fixes 2, 3, 4 independently per chromosome.

    Requires:
        pip install msprime stdpopsim
    """
    try:
        import msprime
        import stdpopsim
    except ImportError:
        raise RuntimeError(
            "msprime and stdpopsim are required. Install with:\n"
            "  pip install msprime stdpopsim"
        )

    vcf_dir   = work_dir / "msprime_vcf"
    plink_dir = work_dir / "plink_chroms"
    vcf_dir.mkdir(exist_ok=True)
    plink_dir.mkdir(exist_ok=True)

    n_parallel = cfg["msprime"].get("n_parallel_chroms", 22)
    log_dir    = cfg["paths"]["logs_dir"]
    chromosomes = cfg["reference"]["chromosomes"]

    # Build argument tuples for each worker
    worker_args = [
        (chrom, cfg, str(work_dir), log_dir)
        for chrom in chromosomes
    ]

    logger.info(
        f"Simulating {len(chromosomes)} chromosomes in parallel "
        f"(n_parallel={n_parallel}) ..."
    )

    # ── FIX 1: Parallel chromosome simulation ──
    with Pool(processes=n_parallel) as pool:
        results = pool.map(_simulate_one_chrom, worker_args)

    logger.info(f"All chromosomes complete: {results}")


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 3: Process Real chrX Data from 1000G — Block Assignment Strategy
#
#  Problem: 1000G has ~2,504 unique chrX samples but we need 75,000.
#  Solution: Assign chrX data in blocks of 2,504, cycling through the real
#  1000G samples repeatedly. Each block is the same real genotype data
#  but relabelled with different SYNTH_ IDs. Sex is preserved — real males
#  always map to male SYNTH_ IDs, real females to female SYNTH_ IDs.
#
#  Result: All 75,000 synthetic individuals have chrX data with consistent
#  SYNTH_ IDs matching sex_assignments_true.txt.
# ─────────────────────────────────────────────────────────────────────────────

def process_chrX(cfg: dict, work_dir: Path, sex_assignments: pd.DataFrame,
                 logger: logging.Logger):
    """
    Extract chrX SNPs from 1000G and assign to ALL synthetic individuals
    using a block assignment strategy:
      - Extract the 2,504 unique 1000G chrX samples once
      - Convert to PLINK .bed with real IDs and sex information
      - Divide 75,000 SYNTH_ IDs into blocks of 2,504
      - For each block: rename the same .bed file with the block's SYNTH_ IDs
      - Merge all blocks into a single chrX .bed with all 75,000 SYNTH_ IDs
    """
    ref_dir      = Path(cfg["paths"]["reference_dir"])
    chrX_vcf     = ref_dir / "ALL.chrX.phase3.vcf.gz"
    plink2       = cfg["paths"]["plink2_bin"]
    plink1       = cfg["paths"].get("plink1_bin", "plink")
    bcftools     = cfg["paths"]["bcftools_bin"]
    tabix_bin    = cfg["paths"].get("tabix_bin", "tabix")
    n_chrX_snps  = cfg["n_snps_chrX"]
    chrX_dir     = work_dir / "chrX"
    chrX_dir.mkdir(exist_ok=True)
    maf_prefilter = cfg["msprime"].get("plink_maf_prefilter", 0.001)

    logger.info("Processing chrX with block assignment strategy ...")

    # ── 1. Get SNP list from chrX VCF ──
    snp_list_file = chrX_dir / "chrX_snp_list.txt"
    cmd = (
        f"{bcftools} view -v snps -m2 -M2 {chrX_vcf} | "
        f"{bcftools} query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > {snp_list_file}"
    )
    run_cmd(cmd, logger)

    snps = pd.read_csv(snp_list_file, sep="\t",
                       names=["CHROM", "POS", "ID", "REF", "ALT"])

    # PAR regions
    PAR1    = (snps["POS"] >= 60001)     & (snps["POS"] <= 2699520)
    PAR2    = (snps["POS"] >= 154931044) & (snps["POS"] <= 155260560)
    non_PAR = snps[~(PAR1 | PAR2)].copy()
    par_snps = snps[PAR1 | PAR2].copy()

    n_non_par = min(int(n_chrX_snps * 0.85), len(non_PAR))
    n_par     = min(n_chrX_snps - n_non_par, len(par_snps))

    sampled_non_par = non_PAR.sample(n=n_non_par, random_state=cfg["sex_distortion"]["seed"])
    sampled_par     = par_snps.sample(n=n_par,    random_state=cfg["sex_distortion"]["seed"])
    sampled_snps    = pd.concat([sampled_non_par, sampled_par]).sort_values("POS")

    # Write SNP position filter file (CHROM:POS) — more reliable than ID
    # because many 1000G chrX variants have ID="." which bcftools cannot match
    snp_id_file  = chrX_dir / "chrX_sampled_snps.txt"
    snp_pos_file = chrX_dir / "chrX_sampled_pos.txt"
    sampled_snps["ID"].to_csv(snp_id_file, index=False, header=False)
    # positions file for bcftools regions filter: CHROM\tPOS (1-based)
    sampled_snps[["CHROM", "POS"]].to_csv(
        snp_pos_file, sep="\t", index=False, header=False
    )
    logger.info(f"chrX: sampled {len(sampled_snps)} SNPs ({n_non_par} non-PAR, {n_par} PAR).")

    # ── 2. Load 1000G sample sex info ──
    sample_info_file = ref_dir / "1000GP_Phase3.sample"
    sample_info = pd.read_csv(sample_info_file, sep=" ", header=0)
    sample_info.columns = sample_info.columns.str.lower()
    real_males   = sample_info[sample_info["sex"].str.lower() == "male"]["id"].tolist()
    real_females = sample_info[sample_info["sex"].str.lower() == "female"]["id"].tolist()
    unique_males   = list(dict.fromkeys(real_males))
    unique_females = list(dict.fromkeys(real_females))
    all_real_ids   = unique_males + unique_females
    block_size     = len(all_real_ids)

    logger.info(
        f"chrX: {len(unique_males)} unique real males, "
        f"{len(unique_females)} unique real females. "
        f"Block size = {block_size}."
    )

    # ── 3. Extract chrX VCF subset (all unique 1000G samples, once) ──
    sample_file = chrX_dir / "chrX_sample_subset.txt"
    pd.Series(all_real_ids).to_csv(sample_file, index=False, header=False)

    chrX_sub_vcf = chrX_dir / "chrX_subset.vcf.gz"
    cmd = (
        f"{bcftools} view "
        f"--samples-file {sample_file} "
        f"--regions-file {snp_pos_file} "   # filter by position, not ID (IDs may be ".")
        f"-v snps -m2 -M2 "                 # biallelic SNPs only
        f"{chrX_vcf} "
        f"-O z -o {chrX_sub_vcf} && "
        f"{tabix_bin} -p vcf {chrX_sub_vcf}"
    )
    run_cmd(cmd, logger)

    # ── 4. Write sex file using real 1000G IDs ──
    sex_file_chrX = chrX_dir / "chrX_sex.txt"
    with open(sex_file_chrX, "w") as sf:
        for real_id in unique_males:
            sf.write(f"{real_id}\t{real_id}\t1\n")
        for real_id in unique_females:
            sf.write(f"{real_id}\t{real_id}\t2\n")

    # ── 5. Convert VCF to base .bed with real 1000G IDs ──
    chrX_base = chrX_dir / "chrX_base"
    cmd = (
        f"{plink2} "
        f"--vcf {chrX_sub_vcf} "
        f"--split-par b37 "
        f"--update-sex {sex_file_chrX} "
        f"--make-bed "
        f"--out {chrX_base} "
        f"--maf {maf_prefilter} "
        f"--max-alleles 2 "
        f"--snps-only "
        f"--set-all-var-ids @:#:$r:$a "
        f"--new-id-max-allele-len 10 missing "
    )
    run_cmd(cmd, logger)

    # ── 5b. Run sex check on base dataset and remove ambiguous samples ──
    # Real 1000G samples with F-statistic in 0.2-0.8 range have ambiguous
    # chrX heterozygosity. When reused 30x in block assignment these propagate
    # to ~25% of synthetic samples causing false sex check failures.
    # We pre-filter to keep only samples with clearly male (F>0.8) or
    # clearly female (F<0.2) chrX profiles.
    chrX_sexcheck_prefix = chrX_dir / "chrX_base_sexcheck"
    cmd = (
        f"{plink1} "
        f"--bfile {chrX_base} "
        f"--check-sex "
        f"--set-hh-missing "   # treat het haploid calls as missing
        f"--allow-extra-chr "
        f"--allow-no-sex "
        f"--out {chrX_sexcheck_prefix} "
    )
    run_cmd(cmd, logger)

    sexcheck_file = Path(str(chrX_sexcheck_prefix) + ".sexcheck")
    ambiguous_ids = set()
    with open(sexcheck_file) as sc:
        next(sc)  # skip header
        for line in sc:
            parts = line.strip().split()
            iid = parts[1]
            try:
                f_stat = float(parts[5])
            except (ValueError, IndexError):
                ambiguous_ids.add(iid)
                continue
            # Keep only samples with clearly male (F>0.8) or female (F<0.2)
            if 0.2 <= f_stat <= 0.8:
                ambiguous_ids.add(iid)

    logger.info(
        f"chrX: {len(ambiguous_ids)} real 1000G samples with ambiguous "
        f"F-statistic (0.2-0.8) will be excluded from block assignment."
    )

    # Remove ambiguous samples from male/female lists
    unique_males   = [x for x in unique_males   if x not in ambiguous_ids]
    unique_females = [x for x in unique_females if x not in ambiguous_ids]

    # Read actual FID values from chrX_base.fam — PLINK2 outputs FID=0
    # so keep file must use the real FID from the .fam, not assume FID=IID
    base_fam_raw = chrX_dir / "chrX_base.fam"
    actual_fid = {}   # IID -> FID as written in .fam
    with open(base_fam_raw) as bf:
        for line in bf:
            parts = line.strip().split()
            actual_fid[parts[1]] = parts[0]  # IID -> FID

    clean_samples = unique_males + unique_females
    clean_keep_file = chrX_dir / "chrX_base_clean_keep.txt"
    with open(clean_keep_file, "w") as kf:
        for iid in clean_samples:
            fid = actual_fid.get(iid, iid)   # fall back to IID if not found
            kf.write(f"{fid}\t{iid}\n")

    logger.info(
        f"chrX: keep file written with {len(clean_samples)} clean samples "
        f"(FID from actual .fam — may be 0 from PLINK2)."
    )

    chrX_base_clean = chrX_dir / "chrX_base_clean"
    cmd = (
        f"{plink1} "
        f"--bfile {chrX_base} "
        f"--keep {clean_keep_file} "
        f"--make-bed "
        f"--out {chrX_base_clean} "
        f"--allow-no-sex "
        f"--allow-extra-chr "
    )
    run_cmd(cmd, logger)

    # Use clean base for all subsequent block operations
    chrX_base = chrX_base_clean
    logger.info(
        f"chrX base after filtering: {len(unique_males)} clean males, "
        f"{len(unique_females)} clean females."
    )

    # Read actual FID/IID from base .fam (FID may be 0 from PLINK2)
    base_fam = Path(str(chrX_base) + ".fam")
    fam_rows = []
    with open(base_fam) as ff:
        for line in ff:
            parts = line.strip().split()
            fam_rows.append((parts[0], parts[1]))  # (FID, IID)

    # Build IID -> sex map from base fam for use in block sex files
    iid_sex_map = {}
    for real_id in unique_males:
        iid_sex_map[real_id] = 1
    for real_id in unique_females:
        iid_sex_map[real_id] = 2

    logger.info(f"chrX base dataset: {len(fam_rows)} clean samples ready for block assignment.")

    # ── 6. Get SYNTH_ IDs split by sex ──
    synth_males   = sex_assignments[sex_assignments["sex"] == 1]["IID"].tolist()
    synth_females = sex_assignments[sex_assignments["sex"] == 2]["IID"].tolist()
    n_synth       = len(sex_assignments)

    # ── 7. Build blocks: assign SYNTH_ IDs sequentially in chunks ──
    # Each block gets the next slice of sex_assignments, split proportionally
    # into males and females matching the real 1000G sex ratio.
    # This guarantees every SYNTH_ ID appears exactly once — no duplicates.
    import math
    n_blocks = math.ceil(n_synth / block_size)
    logger.info(
        f"chrX: assigning {n_synth} synthetic individuals "
        f"across {n_blocks} blocks of size {block_size}."
    )

    blocks_dir = chrX_dir / "blocks"
    blocks_dir.mkdir(exist_ok=True)
    block_prefixes = []

    # Get ordered lists of all SYNTH_ males and females
    all_synth_males   = sex_assignments[sex_assignments["sex"] == 1]["IID"].tolist()
    all_synth_females = sex_assignments[sex_assignments["sex"] == 2]["IID"].tolist()
    male_ptr   = 0   # pointer into all_synth_males
    female_ptr = 0   # pointer into all_synth_females

    for b in range(n_blocks):
        # How many SYNTH_ IDs remain for this block
        synth_remaining = n_synth - b * block_size
        this_block_size = min(block_size, synth_remaining)

        # Split proportionally by sex ratio from real 1000G data
        n_block_males   = round(this_block_size * len(unique_males)   / block_size)
        n_block_females = this_block_size - n_block_males

        # Take the next slice of SYNTH_ IDs for each sex
        block_synth_males   = all_synth_males  [male_ptr   : male_ptr   + n_block_males]
        block_synth_females = all_synth_females[female_ptr : female_ptr + n_block_females]
        male_ptr   += n_block_males
        female_ptr += n_block_females

        if not block_synth_males and not block_synth_females:
            break

        # Build rename file: (actual_FID, real_IID) -> (SYNTH_ID, SYNTH_ID)
        # Real males map to SYNTH_ males, real females to SYNTH_ females
        iid_to_synth = {}
        for real_id, synth_id in zip(unique_males,   block_synth_males):
            iid_to_synth[real_id] = synth_id
        for real_id, synth_id in zip(unique_females, block_synth_females):
            iid_to_synth[real_id] = synth_id

        rename_file = blocks_dir / f"block_{b:04d}_rename.txt"
        keep_file   = blocks_dir / f"block_{b:04d}_keep.txt"
        n_renamed = 0
        with open(rename_file, "w") as rf, open(keep_file, "w") as kf:
            for fid, iid in fam_rows:
                if iid in iid_to_synth:
                    synth_id = iid_to_synth[iid]
                    rf.write(f"{fid}\t{iid}\t{synth_id}\t{synth_id}\n")
                    # keep file uses NEW SYNTH_ IDs (after rename)
                    kf.write(f"{synth_id}\t{synth_id}\n")
                    n_renamed += 1

        # Step A: apply rename — gives all 2504 real samples new SYNTH_ IDs
        # (samples not in rename_file keep their real HG*/NA* IDs)
        block_renamed = blocks_dir / f"block_{b:04d}_renamed"
        cmd = (
            f"{plink1} "
            f"--bfile {chrX_base} "
            f"--update-ids {rename_file} "
            f"--make-bed "
            f"--out {block_renamed} "
            f"--allow-no-sex "
            f"--allow-extra-chr "
        )
        run_cmd(cmd, logger)

        # Step B: keep ONLY the renamed SYNTH_ samples
        # This drops the un-renamed real IDs (HG*/NA*) from the block
        # so they do not appear in the final merged dataset
        block_out = blocks_dir / f"block_{b:04d}"
        cmd = (
            f"{plink1} "
            f"--bfile {block_renamed} "
            f"--keep {keep_file} "
            f"--make-bed "
            f"--out {block_out} "
            f"--allow-no-sex "
            f"--allow-extra-chr "
        )
        run_cmd(cmd, logger)
        block_prefixes.append(str(block_out))

        logger.info(
            f"chrX block {b+1}/{n_blocks}: "
            f"{n_renamed}/{len(fam_rows)} samples renamed and kept."
        )

    logger.info(f"chrX: {len(block_prefixes)} blocks created. Merging ...")

    # ── 8. Combine blocks into single chrX dataset ──
    # PLINK --merge stacks VARIANTS not SAMPLES when all blocks share
    # the same variant IDs — giving N_blocks x N_variants instead of
    # N_total_samples x N_variants. We must combine blocks directly.
    #
    # Strategy: load each block's .bed into a numpy array of shape
    # (N_variants, N_samples_in_block), concatenate along axis=1
    # (samples), then write a single combined .bed.
    #
    # PLINK .bed SNP-major format:
    #   3-byte magic header
    #   For each variant: ceil(N/4) bytes encoding N 2-bit genotypes
    #   Bit order within each byte: LSB first (samples 0,1,2,3 in bits 0-1,2-3,4-5,6-7)
    chrX_plink = chrX_dir / "chrX"
    import shutil

    if len(block_prefixes) == 1:
        for ext in [".bed", ".bim", ".fam"]:
            shutil.copy(str(block_prefixes[0]) + ext, str(chrX_plink) + ext)
        logger.info("chrX: single block, copied directly.")
    else:
        logger.info(
            f"chrX: combining {len(block_prefixes)} blocks "
            f"using numpy-based .bed concatenation ..."
        )

        PLINK_MAGIC = b"\x6c\x1b\x01"

        # Count samples per block and total variants (same for all blocks)
        n_variants = sum(1 for _ in open(str(block_prefixes[0]) + ".bim"))
        n_per_block = [
            sum(1 for _ in open(str(p) + ".fam"))
            for p in block_prefixes
        ]
        n_total = sum(n_per_block)
        logger.info(
            f"chrX: {n_variants} variants x {n_total} total samples "
            f"({len(block_prefixes)} blocks)."
        )

        # Read each block .bed into a 2D numpy array: shape (N_variants, N_samples)
        # PLINK stores 4 genotypes per byte; unpackbits gives us 8 bits per byte.
        # We unpack then reshape to (N_variants, N_bytes_per_row * 4),
        # then slice to exact N_samples columns.
        block_arrays = []
        for prefix, n_samp in zip(block_prefixes, n_per_block):
            bed_path = str(prefix) + ".bed"
            with open(bed_path, "rb") as fh:
                magic = fh.read(3)
                if magic != PLINK_MAGIC:
                    raise RuntimeError(
                        f"Bad .bed magic in {bed_path}: {magic!r}"
                    )
                raw = np.frombuffer(fh.read(), dtype=np.uint8)

            n_bytes_row = math.ceil(n_samp / 4)
            # raw has shape (N_variants * n_bytes_row,)
            raw = raw.reshape(n_variants, n_bytes_row)

            # Unpack 8 bits per byte → shape (N_variants, n_bytes_row * 8)
            bits = np.unpackbits(raw, axis=1, bitorder="little")

            # Each genotype is 2 bits; reshape to (N_variants, n_bytes_row*4, 2)
            genos = bits[:, : n_bytes_row * 4 * 2].reshape(
                n_variants, n_bytes_row * 4, 2
            )
            # Trim padding columns to exact sample count
            genos = genos[:, :n_samp, :]  # (N_variants, n_samp, 2)
            block_arrays.append(genos)

        # Concatenate all blocks along the samples axis
        # Result shape: (N_variants, N_total, 2)
        combined = np.concatenate(block_arrays, axis=1)

        # Re-pack into bytes
        n_bytes_out = math.ceil(n_total / 4)
        # Pad to multiple of 4 samples with zeros — PLINK ignores padding bits
        pad = n_bytes_out * 4 - n_total
        if pad > 0:
            padding = np.zeros((n_variants, pad, 2), dtype=np.uint8)
            combined = np.concatenate([combined, padding], axis=1)

        # Flatten 2-bit pairs back to bits: (N_variants, n_bytes_out*8)
        flat_bits = combined.reshape(n_variants, n_bytes_out * 8)
        packed = np.packbits(flat_bits.astype(np.uint8), axis=1, bitorder="little")
        # packed shape: (N_variants, n_bytes_out)

        chrX_bed = Path(str(chrX_plink) + ".bed")
        with open(chrX_bed, "wb") as out:
            out.write(PLINK_MAGIC)
            out.write(packed.tobytes())

        # Stack .fam files
        chrX_fam = Path(str(chrX_plink) + ".fam")
        with open(chrX_fam, "w") as out:
            for prefix in block_prefixes:
                with open(str(prefix) + ".fam") as fh:
                    out.write(fh.read())

        # .bim is identical across all blocks
        shutil.copy(str(block_prefixes[0]) + ".bim", str(chrX_plink) + ".bim")

        logger.info(
            f"chrX: numpy concatenation complete — "
            f"{n_total} samples x {n_variants} variants."
        )

    # ── 9. Verify final chrX fam ──
    final_fam = chrX_dir / "chrX.fam"
    n_final = sum(1 for _ in open(final_fam))
    with open(final_fam) as ff:
        first_iid = ff.readline().strip().split()[1]

    logger.info(
        f"chrX processing complete: {n_final} samples in final dataset. "
        f"First IID: {first_iid}"
    )
    if not first_iid.startswith("SYNTH_"):
        logger.warning(
            f"chrX: First IID is {first_iid!r} — expected SYNTH_ format. "
            f"Check block rename files in {blocks_dir}."
        )

    return chrX_plink


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 4: Assign Sex + Merge All Chromosomes
# ─────────────────────────────────────────────────────────────────────────────

def assign_sex(cfg: dict, work_dir: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Create sex assignments for N individuals based on configured female_fraction.
    Returns a DataFrame with columns: FID, IID, sex (1=male, 2=female)
    """
    n        = cfg["n_individuals"]
    n_female = int(n * cfg["sex_ratio"]["female_fraction"])
    n_male   = n - n_female

    rng        = np.random.default_rng(cfg["sex_distortion"]["seed"])
    sex_labels = np.array([2] * n_female + [1] * n_male)
    rng.shuffle(sex_labels)

    iids   = [f"SYNTH_{i+1:07d}" for i in range(n)]
    sex_df = pd.DataFrame({"FID": iids, "IID": iids, "sex": sex_labels})

    sex_file = work_dir / "sex_assignments_true.txt"
    sex_df.to_csv(sex_file, sep="\t", index=False)
    logger.info(f"Sex assigned: {n_female} females, {n_male} males.")
    return sex_df


def _fix_fid_in_fam(out_prefix: Path, plink1: str, logger: logging.Logger):
    """
    Fix FID=0 issue after PLINK1.9 merge.
    PLINK1.9 sets FID=0 for all samples when merging files without family
    structure. This reads the merged .fam, sets FID=IID for every sample,
    writes an --update-ids file, and applies it with plink --update-ids.
    """
    fam_file = Path(str(out_prefix) + ".fam")
    if not fam_file.exists():
        logger.warning(f"FAM file not found for FID fix: {fam_file}")
        return

    # Read current fam — columns: FID IID PAT MAT SEX PHENO
    rows = []
    with open(fam_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                rows.append(parts)

    # Build update-ids file: OLD_FID OLD_IID NEW_FID NEW_IID
    # Set NEW_FID = IID so FID matches IID throughout
    update_ids_file = Path(str(out_prefix) + "_update_fid.txt")
    needs_fix = any(r[0] != r[1] for r in rows)

    if not needs_fix:
        logger.info("FID already matches IID — no fix needed.")
        return

    with open(update_ids_file, "w") as uf:
        for r in rows:
            old_fid, iid = r[0], r[1]
            uf.write(f"{old_fid}\t{iid}\t{iid}\t{iid}\n")

    # Apply the fix using plink --update-ids in-place
    fixed_prefix = Path(str(out_prefix) + "_fidfixed")
    cmd = (
        f"{plink1} "
        f"--bfile {out_prefix} "
        f"--update-ids {update_ids_file} "
        f"--make-bed "
        f"--out {fixed_prefix} "
        f"--allow-no-sex "
        f"--allow-extra-chr "
    )
    import subprocess as _sp
    result = _sp.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        logger.warning(f"FID fix failed: {result.stderr}. Keeping FID=0.")
        return

    # Replace original files with fixed files
    for ext in [".bed", ".bim", ".fam"]:
        fixed_file   = Path(str(fixed_prefix) + ext)
        original_file = Path(str(out_prefix) + ext)
        if fixed_file.exists():
            fixed_file.replace(original_file)

    # Clean up temp files
    update_ids_file.unlink(missing_ok=True)
    for ext in [".log", ".nosex"]:
        tmp = Path(str(fixed_prefix) + ext)
        if tmp.exists():
            tmp.unlink()

    logger.info(f"FID fixed: set FID=IID for {len(rows)} samples.")


def merge_all_chromosomes(cfg: dict, work_dir: Path, logger: logging.Logger) -> Path:
    """
    Merge per-chromosome PLINK .bed/.bim/.fam files into a single dataset.
    Uses PLINK1.9 (plink) --merge-list which natively supports .bed merging.
    PLINK2 --pmerge-list requires .pgen format so cannot be used here.
    """
    # Use plink1.9 for merging — plink2 does not support --merge-list for .bed files
    plink1     = cfg["paths"].get("plink1_bin", "plink")
    plink_dir  = work_dir / "plink_chroms"
    merged_dir = work_dir / "merged"
    merged_dir.mkdir(exist_ok=True)

    # Build merge list — plink1.9 --merge-list expects one prefix per line
    merge_list = merged_dir / "merge_list.txt"
    all_prefixes = []
    for chrom in cfg["reference"]["chromosomes"]:
        pfile = plink_dir / f"chr{chrom}"
        if (plink_dir / f"chr{chrom}.bed").exists():
            all_prefixes.append(str(pfile))
    chrX_pfile = work_dir / "chrX" / "chrX"
    if (chrX_pfile.parent / "chrX.bed").exists():
        all_prefixes.append(str(chrX_pfile))

    if len(all_prefixes) == 0:
        raise RuntimeError("No .bed files found to merge. Run --step simulate and --step chrX first.")

    logger.info(f"Merging {len(all_prefixes)} chromosome files ...")

    # plink1.9 --merge-list: first file is the primary, rest go in the list file
    primary = all_prefixes[0]
    with open(merge_list, "w") as f:
        for prefix in all_prefixes[1:]:
            f.write(prefix + "\n")

    out_prefix   = merged_dir / "synthetic_full"
    missnp_file  = Path(str(out_prefix) + "-merge.missnp")
    exclude_file = merged_dir / "all_excluded_multiallelic.txt"

    # ── Pass 1: Probe merge to discover all multiallelic variants ──
    # PLINK1.9 --exclude during merge only filters the PRIMARY file,
    # NOT the files in the merge list. So we must pre-filter EACH
    # chromosome file individually before the final merge.
    logger.info("Pass 1: probe merge to discover multiallelic variants ...")
    cmd = (
        f"{plink1} "
        f"--bfile {primary} "
        f"--merge-list {merge_list} "
        f"--make-bed "
        f"--out {out_prefix} "
        f"--allow-no-sex "
        f"--allow-extra-chr "
    )
    result = run_cmd(cmd, logger, check=False)

    # ── If merge succeeded with no conflicts, we are done ──
    if result.returncode == 0 and (
        not missnp_file.exists() or missnp_file.stat().st_size == 0
    ):
        logger.info("Merge succeeded on first attempt with no multiallelic conflicts.")
        _fix_fid_in_fam(out_prefix, plink1, logger)
        logger.info(f"Merged dataset: {out_prefix}")
        return out_prefix

    # ── Collect all multiallelic variant IDs from missnp file ──
    if not missnp_file.exists() or missnp_file.stat().st_size == 0:
        raise RuntimeError(
            f"Merge failed with no missnp file:\n{result.stderr}"
        )

    excluded_snps = set()
    with open(missnp_file) as mf:
        for line in mf:
            snp = line.strip()
            if snp:
                excluded_snps.add(snp)

    logger.warning(
        f"Found {len(excluded_snps)} multiallelic variants. "
        f"Pre-filtering ALL chromosome files before re-merging ..."
    )

    # Write the exclusion list
    with open(exclude_file, "w") as ef:
        for snp in sorted(excluded_snps):
            ef.write(snp + "\n")

    # ── Pass 2: Pre-filter EVERY chromosome file individually ──
    # This is the correct fix — PLINK1.9 --exclude during merge only
    # applies to the primary file. We must filter each file separately
    # so the multiallelic variants are gone from ALL inputs before merging.
    filtered_dir = merged_dir / "filtered_chroms"
    filtered_dir.mkdir(exist_ok=True)
    filtered_prefixes = []

    for prefix in all_prefixes:
        chrom_name   = Path(prefix).name
        filtered_out = filtered_dir / chrom_name
        cmd = (
            f"{plink1} "
            f"--bfile {prefix} "
            f"--exclude {exclude_file} "
            f"--make-bed "
            f"--out {filtered_out} "
            f"--allow-no-sex "
            f"--allow-extra-chr "
        )
        run_cmd(cmd, logger)
        filtered_prefixes.append(str(filtered_out))

    logger.info(
        f"Pre-filtered {len(filtered_prefixes)} chromosome files. "
        f"Excluded {len(excluded_snps)} multiallelic variants from each."
    )

    # ── Pass 3: Final merge using pre-filtered files ──
    filtered_primary   = filtered_prefixes[0]
    filtered_merge_list = merged_dir / "filtered_merge_list.txt"
    with open(filtered_merge_list, "w") as f:
        for prefix in filtered_prefixes[1:]:
            f.write(prefix + "\n")

    cmd = (
        f"{plink1} "
        f"--bfile {filtered_primary} "
        f"--merge-list {filtered_merge_list} "
        f"--make-bed "
        f"--out {out_prefix} "
        f"--allow-no-sex "
        f"--allow-extra-chr "
    )
    run_cmd(cmd, logger)

    logger.info(
        f"Merge complete. Excluded {len(excluded_snps)} multiallelic variants total."
    )
    logger.info(f"Excluded SNPs list: {exclude_file}")

    # ── Fix FID=0 issue ──
    # PLINK1.9 sets FID=0 for all samples when merging files with no
    # family structure. We fix this by setting FID=IID for all samples
    # so the merged .fam is consistent with sex_assignments_true.txt.
    _fix_fid_in_fam(out_prefix, plink1, logger)

    logger.info(f"Merged dataset: {out_prefix}")
    return out_prefix


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 5: Simulate Phenotypes (calls R/PhenotypeSimulator)
# ─────────────────────────────────────────────────────────────────────────────

def simulate_phenotypes(cfg: dict, merged_prefix: Path, work_dir: Path,
                        logger: logging.Logger) -> Path:
    """
    Call the R phenotype simulation script (simulate_phenotypes.R).
    Outputs phenotypes.tsv, causal_snps.tsv, covariates.tsv.
    """
    pheno_dir      = work_dir / "phenotypes"
    pheno_dir.mkdir(exist_ok=True)
    r_script       = Path(__file__).parent / "simulate_phenotypes.R"
    pheno_cfg_json = pheno_dir / "pheno_config.json"

    pheno_params = {
        "plink_prefix":   str(merged_prefix),
        "out_dir":        str(pheno_dir),
        "n_causal_snps":  cfg["phenotype"]["n_causal_snps"],
        "heritability":   cfg["phenotype"]["heritability"],
        "pheno_type":     cfg["phenotype"]["type"],
        "pheno_subtype":  cfg["phenotype"].get("pheno_subtype", "FEV1"),
        "prevalence":     cfg["phenotype"]["prevalence"],
        "strat_variance": cfg["phenotype"]["stratification_variance"],
        "seed":           cfg["sex_distortion"]["seed"],
    }
    with open(pheno_cfg_json, "w") as f:
        json.dump(pheno_params, f, indent=2)

    cmd = f"Rscript {r_script} --config {pheno_cfg_json}"
    run_cmd(cmd, logger)
    logger.info("Phenotype simulation complete.")
    return pheno_dir


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 5b: Simulate Realistic Data Quality
#  Introduces heterogeneous missingness into synthetic genotype data so that
#  QC steps (--mind, --geno, heterozygosity filtering) actually remove samples
#  and variants, as in real GWAS data.
#
#  Three tiers of missingness are simulated:
#    - Background: low uniform missingness across all samples/SNPs
#    - Poor quality samples: a small fraction with high missingness (>5%)
#    - Poor quality SNPs: a small fraction with high missingness (>5%)
#
#  Missingness is introduced by directly editing the .bed binary file —
#  randomly setting genotype calls to 0b01 (PLINK missing code).
# ─────────────────────────────────────────────────────────────────────────────

def simulate_data_quality(
    plink_prefix: Path,
    cfg: dict,
    dataset_id: int,
    logger: logging.Logger,
) -> None:
    """
    Introduce realistic heterogeneous missingness into a PLINK .bed file.

    Missingness model (all rates configurable in config.yaml):
      background_missing_rate : applied uniformly to all genotypes
      poor_sample_fraction    : fraction of samples flagged as poor quality
      poor_sample_miss_rate   : missingness rate for poor-quality samples
      poor_snp_fraction       : fraction of SNPs flagged as poor quality
      poor_snp_miss_rate      : missingness rate for poor-quality SNPs

    The function edits the .bed file in-place using numpy.
    Ground-truth lists of poor-quality samples/SNPs are written alongside
    the replicate for benchmarking purposes.
    """
    qc_cfg  = cfg.get("data_quality", {})
    if not qc_cfg.get("enabled", True):
        logger.info("Data quality simulation disabled in config — skipping.")
        return

    bg_rate         = qc_cfg.get("background_missing_rate", 0.005)
    poor_samp_frac  = qc_cfg.get("poor_sample_fraction",    0.03)
    poor_samp_rate  = qc_cfg.get("poor_sample_miss_rate",   0.15)
    poor_snp_frac   = qc_cfg.get("poor_snp_fraction",       0.02)
    poor_snp_rate   = qc_cfg.get("poor_snp_miss_rate",      0.10)
    seed            = cfg["sex_distortion"]["seed"] + dataset_id * 9999

    rng = np.random.default_rng(seed)

    bed_path = Path(str(plink_prefix) + ".bed")
    bim_path = Path(str(plink_prefix) + ".bim")
    fam_path = Path(str(plink_prefix) + ".fam")

    # Read dimensions
    n_variants = sum(1 for _ in open(bim_path))
    n_samples  = sum(1 for _ in open(fam_path))
    n_bytes_row = math.ceil(n_samples / 4)

    logger.info(
        f"Dataset {dataset_id:03d}: simulating data quality "
        f"({n_variants} variants x {n_samples} samples) ..."
    )

    # Read full .bed into memory
    PLINK_MAGIC = b"\x6c\x1b\x01"
    with open(bed_path, "rb") as f:
        magic = f.read(3)
        if magic != PLINK_MAGIC:
            raise RuntimeError(f"Bad .bed magic: {magic!r}")
        raw = np.frombuffer(f.read(), dtype=np.uint8).copy()

    # Reshape to (N_variants, N_bytes_row)
    raw = raw.reshape(n_variants, n_bytes_row)

    # Unpack to (N_variants, N_samples) 2-bit genotype array
    bits  = np.unpackbits(raw, axis=1, bitorder="little")
    genos = bits[:, : n_bytes_row * 4 * 2].reshape(
        n_variants, n_bytes_row * 4, 2
    )[:, :n_samples, :]   # (N_variants, N_samples, 2)

    # PLINK missing genotype = 0b01
    # PLINK missing code = 0b01, i.e. bit0=1, bit1=0 in LSB-first encoding
    # numpy unpackbits(bitorder="little") gives (bit0, bit1) as [bit0_val, bit1_val]
    # So missing = [1, 0], NOT [0, 1] which would be heterozygous (0b10)
    MISSING = np.array([1, 0], dtype=np.uint8)

    # ── Select poor-quality samples ──
    n_poor_samp = max(1, int(n_samples * poor_samp_frac))
    poor_samp_idx = rng.choice(n_samples, size=n_poor_samp, replace=False)

    # ── Select poor-quality SNPs ──
    n_poor_snp = max(1, int(n_variants * poor_snp_frac))
    poor_snp_idx = rng.choice(n_variants, size=n_poor_snp, replace=False)

    # ── Apply background missingness to ALL genotypes ──
    miss_mask = rng.random((n_variants, n_samples)) < bg_rate
    genos[miss_mask] = MISSING

    # ── Apply elevated missingness to poor-quality samples ──
    for si in poor_samp_idx:
        miss_mask_s = rng.random(n_variants) < poor_samp_rate
        genos[miss_mask_s, si] = MISSING

    # ── Apply elevated missingness to poor-quality SNPs ──
    for vi in poor_snp_idx:
        miss_mask_v = rng.random(n_samples) < poor_snp_rate
        genos[vi, miss_mask_v] = MISSING

    # Repack back to .bed bytes
    pad = n_bytes_row * 4 - n_samples
    if pad > 0:
        padding = np.zeros((n_variants, pad, 2), dtype=np.uint8)
        genos_padded = np.concatenate([genos, padding], axis=1)
    else:
        genos_padded = genos

    flat_bits = genos_padded.reshape(n_variants, n_bytes_row * 8)
    packed    = np.packbits(
        flat_bits.astype(np.uint8), axis=1, bitorder="little"
    )

    # Write back
    with open(bed_path, "wb") as f:
        f.write(PLINK_MAGIC)
        f.write(packed.tobytes())

    # ── Write ground-truth quality files for benchmarking ──
    rep_dir = plink_prefix.parent

    # Poor quality samples
    sample_iids = []
    with open(fam_path) as f:
        for line in f:
            sample_iids.append(line.strip().split()[1])
    poor_samp_file = rep_dir / "ground_truth_poor_samples.txt"
    with open(poor_samp_file, "w") as f:
        for idx in sorted(poor_samp_idx):
            iid = sample_iids[idx]
            f.write(f"{iid}\t{iid}\n")

    # Poor quality SNPs
    snp_ids = []
    with open(bim_path) as f:
        for line in f:
            snp_ids.append(line.strip().split()[1])
    poor_snp_file = rep_dir / "ground_truth_poor_snps.txt"
    with open(poor_snp_file, "w") as f:
        for idx in sorted(poor_snp_idx):
            f.write(snp_ids[idx] + "\n")

    logger.info(
        f"Dataset {dataset_id:03d}: quality simulation complete. "
        f"{n_poor_samp} poor-quality samples ({poor_samp_frac*100:.1f}%, "
        f"miss_rate={poor_samp_rate*100:.0f}%), "
        f"{n_poor_snp} poor-quality SNPs ({poor_snp_frac*100:.1f}%, "
        f"miss_rate={poor_snp_rate*100:.0f}%), "
        f"background_rate={bg_rate*100:.1f}%."
    )
    logger.info(f"  Ground-truth poor samples -> {poor_samp_file}")
    logger.info(f"  Ground-truth poor SNPs    -> {poor_snp_file}")


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 6: Introduce Sex Distortion (core of the benchmarking design)
# ─────────────────────────────────────────────────────────────────────────────

def introduce_sex_distortion(cfg: dict, merged_prefix: Path, work_dir: Path,
                              pheno_dir: Path, sex_df: pd.DataFrame,
                              dataset_id: int, logger: logging.Logger) -> dict:
    """
    For a given dataset_id, draw a random distortion percentage from Uniform(0%, 5%),
    flip sex labels for that fraction of individuals, and write a complete
    self-contained PLINK dataset for that replicate.
    """
    plink2   = cfg["paths"]["plink2_bin"]
    dist_cfg = cfg["sex_distortion"]

    rng = np.random.default_rng(dist_cfg["seed"] + dataset_id * 1000)
    if dist_cfg["mode"] == "random":
        distortion_pct = float(rng.uniform(dist_cfg["min_pct"], dist_cfg["max_pct"]))
    else:
        distortion_pct = float(dist_cfg["fixed_pct"])

    n_to_distort = int(len(sex_df) * distortion_pct / 100)
    logger.info(
        f"Dataset {dataset_id:03d}: distortion = {distortion_pct:.4f}% "
        f"({n_to_distort} individuals)"
    )

    distorted_sex = sex_df[["FID", "IID", "sex"]].copy()
    if n_to_distort > 0:
        flip_indices = rng.choice(len(distorted_sex), size=n_to_distort, replace=False)
        distorted_sex.loc[flip_indices, "sex"] = distorted_sex.loc[
            flip_indices, "sex"
        ].map({1: 2, 2: 1})

    rep_dir = work_dir / "replicates" / f"dataset_{dataset_id:03d}"
    rep_dir.mkdir(parents=True, exist_ok=True)

    rep_prefix = rep_dir / "genotypes"
    plink1 = cfg["paths"].get("plink1_bin", "plink")

    # ── Copy merged dataset to replicate directory ──
    import shutil
    for ext in [".bed", ".bim", ".fam"]:
        shutil.copy(str(merged_prefix) + ext, str(rep_prefix) + ext)

    # ── Apply distorted sex directly to .fam file in Python ──
    # This bypasses PLINK --update-sex entirely, avoiding all FID/IID
    # matching issues. We build a lookup IID -> distorted_sex and
    # rewrite column 5 of the .fam file directly.
    iid_to_sex = dict(zip(
        distorted_sex["IID"].tolist(),
        distorted_sex["sex"].tolist()
    ))

    rep_fam = Path(str(rep_prefix) + ".fam")
    fam_lines = rep_fam.read_text().splitlines()
    updated_lines = []
    n_updated = 0
    for line in fam_lines:
        if not line.strip():
            updated_lines.append(line)
            continue
        parts = line.split()
        iid = parts[1]
        if iid in iid_to_sex:
            parts[4] = str(iid_to_sex[iid])
            n_updated += 1
        updated_lines.append("\t".join(parts))
    rep_fam.write_text("\n".join(updated_lines) + "\n")
    logger.info(
        f"Dataset {dataset_id:03d}: sex written directly to .fam "
        f"for {n_updated}/{len(fam_lines)} samples."
    )

    # ── Verify sex was actually updated ──
    rep_fam = Path(str(rep_prefix) + ".fam")
    sex_counts = {}
    with open(rep_fam) as rf:
        for line in rf:
            s = line.strip().split()[4]
            sex_counts[s] = sex_counts.get(s, 0) + 1
    n_unknown = sex_counts.get("0", 0)
    n_male    = sex_counts.get("1", 0)
    n_female  = sex_counts.get("2", 0)
    logger.info(
        f"Dataset {dataset_id:03d} sex check: "
        f"{n_male} males, {n_female} females, {n_unknown} unknown."
    )
    if n_unknown > 0:
        logger.warning(
            f"Dataset {dataset_id:03d}: {n_unknown} samples have unknown sex (col5=0). "
            f"--update-sex may have failed. Check {sex_update_file}."
        )

    # ── Simulate realistic data quality (missingness) ──
    simulate_data_quality(rep_prefix, cfg, dataset_id, logger)

    for fname in ["phenotypes.tsv", "causal_snps.tsv", "covariates.tsv"]:
        src = pheno_dir / fname
        dst = rep_dir / fname
        if src.exists() and not dst.exists():
            os.symlink(src.resolve(), dst)

    distorted_iids = distorted_sex.loc[
        flip_indices if n_to_distort > 0 else [], "IID"
    ].tolist() if n_to_distort > 0 else []

    metadata = {
        "dataset_id":     dataset_id,
        "distortion_pct": round(distortion_pct, 6),
        "n_distorted":    n_to_distort,
        "distorted_iids": distorted_iids,
        "seed_used":      int(dist_cfg["seed"] + dataset_id * 1000),
        "plink_prefix":   str(rep_prefix),
        "phenotype_dir":  str(pheno_dir),
    }
    meta_file = rep_dir / "metadata.json"
    with open(meta_file, "w") as f:
        json.dump(metadata, f, indent=2)

    logger.info(f"Dataset {dataset_id:03d} written -> {rep_dir}")
    return metadata


# ─────────────────────────────────────────────────────────────────────────────
#  STEP 7: Summary Manifest
# ─────────────────────────────────────────────────────────────────────────────

def write_manifest(all_metadata: list, work_dir: Path, logger: logging.Logger):
    """Write a summary TSV of all generated replicates."""
    manifest = pd.DataFrame([
        {
            "dataset_id":     m["dataset_id"],
            "distortion_pct": m["distortion_pct"],
            "n_distorted":    m["n_distorted"],
            "plink_prefix":   m["plink_prefix"],
        }
        for m in all_metadata
    ])
    out = work_dir / "replicates_manifest.tsv"
    manifest.to_csv(out, sep="\t", index=False)
    logger.info(f"Manifest written -> {out}")
    logger.info("\n" + manifest.to_string(index=False))


# ─────────────────────────────────────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="GWAS Synthetic Data Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Single run with default seed from config\n"
            "  python gwas_pipeline_master.py --config config.yaml --step all\n"
            "\n"
            "  # Three independent datasets with different seeds\n"
            "  python gwas_pipeline_master.py --config config.yaml --seed 42 --step all\n"
            "  python gwas_pipeline_master.py --config config.yaml --seed 123 --step all\n"
            "  python gwas_pipeline_master.py --config config.yaml --seed 456 --step all\n"
        )
    )
    parser.add_argument("--config",     default="config.yaml",
                        help="Path to config.yaml")
    parser.add_argument("--step",       default="all",
                        choices=["all", "download", "simulate", "chrX",
                                 "merge", "phenotype", "distort"],
                        help="Run a specific step only")
    parser.add_argument("--dataset-id", type=int, default=None,
                        help="Run distortion for a single dataset ID (for parallelisation)")
    parser.add_argument("--seed",       type=int, default=None,
                        help=(
                            "Override master seed from config.yaml. "
                            "Each unique seed produces a fully independent "
                            "synthetic population with different genotypes and phenotypes. "
                            "Output is saved to a seed-specific subdirectory e.g. "
                            "gwas_synthetic_output/seed_42/"
                        ))
    args = parser.parse_args()

    cfg = load_config(args.config)

    # ── Seed override ──
    # If --seed is provided it overrides the seed in config.yaml and
    # redirects all output to a seed-specific subdirectory so multiple
    # independent runs do not overwrite each other.
    if args.seed is not None:
        cfg["sex_distortion"]["seed"] = args.seed
        base_dir = Path(cfg["paths"]["working_dir"])
        seed_dir = base_dir / f"seed_{args.seed}"
        cfg["paths"]["working_dir"] = str(seed_dir)
        cfg["paths"]["logs_dir"]    = str(seed_dir / "logs")

    work_dir = Path(cfg["paths"]["working_dir"])
    work_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging(cfg["paths"]["logs_dir"], args.dataset_id)
    logger.info("=" * 70)
    logger.info("  GWAS Synthetic Data Pipeline - Starting")
    logger.info(f"  Config:     {args.config}")
    logger.info(f"  Step:       {args.step}")
    logger.info(f"  Seed:       {cfg['sex_distortion']['seed']}")
    logger.info(f"  Output dir: {work_dir}")
    logger.info("=" * 70)

    run_all = (args.step == "all")

    # Step 1: Download reference (shared across all seeds — never re-downloaded)
    if run_all or args.step == "download":
        download_1000g_reference(cfg, logger)

    # Step 2: Simulate autosomes (msprime, parallel, 50MB regions, thinned VCF)
    if run_all or args.step == "simulate":
        simulate_autosomes_hapgen2(cfg, work_dir, logger)

    # Step 3+4: Sex assignment + chrX processing
    sex_df = assign_sex(cfg, work_dir, logger)
    if run_all or args.step == "chrX":
        process_chrX(cfg, work_dir, sex_df, logger)

    # Step 5: Merge chromosomes
    if run_all or args.step == "merge":
        merged_prefix = merge_all_chromosomes(cfg, work_dir, logger)
    else:
        merged_prefix = work_dir / "merged" / "synthetic_full"

    # Step 6: Simulate phenotypes
    if run_all or args.step == "phenotype":
        pheno_dir = simulate_phenotypes(cfg, merged_prefix, work_dir, logger)
    else:
        pheno_dir = work_dir / "phenotypes"

    # Step 7: Generate replicates with sex distortion
    if run_all or args.step == "distort":
        n_datasets   = cfg["sex_distortion"]["n_datasets"]
        all_metadata = []

        if args.dataset_id is not None:
            meta = introduce_sex_distortion(
                cfg, merged_prefix, work_dir, pheno_dir, sex_df,
                args.dataset_id, logger
            )
            all_metadata.append(meta)
        else:
            for ds_id in range(1, n_datasets + 1):
                meta = introduce_sex_distortion(
                    cfg, merged_prefix, work_dir, pheno_dir, sex_df,
                    ds_id, logger
                )
                all_metadata.append(meta)

        write_manifest(all_metadata, work_dir, logger)

    logger.info("=" * 70)
    logger.info(f"  Pipeline complete — seed={cfg['sex_distortion']['seed']}")
    logger.info(f"  Output: {work_dir}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
