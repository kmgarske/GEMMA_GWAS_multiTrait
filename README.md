# Use the GEMMA software to run GWAS across multiple traits

A Snakemake workflow for running a GEMMA genome-wide association study (GWAS) across **many traits at once** using a linear mixed model with **leave-one-chromosome-out (LOCO)** kinship matrices.

For every trait in an input trait matrix, the pipeline fits a univariate LMM (`gemma -lmm 4`) on each autosome (chr1–chr22), correcting for relatedness/population structure with a LOCO kinship matrix and for user-supplied covariates plus genotype principal components.

## Overview

This pipeline:

- Extracts a per-trait phenotype vector from a multi-trait matrix and matches it to the genotyped samples
- Filters genotypes by MAF, missingness, and call rate, and LD-prunes a copy for kinship/PCA
- Computes genotype principal components for population structure correction
- Builds **LOCO kinship matrices** (one per chromosome, chr1–chr22) for the mixed model
- Assembles a covariate matrix combining genotype PCs and user covariates
- Runs a per-trait, per-chromosome GEMMA LMM and writes association results

## Installation

### Dependencies

Create a conda environment from the provided spec:

```bash
conda env create -f workflow/envs/allEnv.yml
conda activate <env_name>
```

Key dependencies (see `workflow/envs/allEnv.yml`):

- snakemake = 9.6.0 (plus `snakemake-executor-plugin-slurm`)
- gemma = 0.98.3
- plink = 1.90b6.21
- bcftools = 1.20
- bedtools = 2.31.1
- gatk4 = 4.5.0.0
- R = 4.3.3 with tidyverse 2.0.0
- pandas = 2.1.4

## Input File Formats

### 1. Genotype files (`inputBfile`)

PLINK binary format (`.bed`, `.bim`, `.fam`). Sample IDs in the `.fam` file should be the unique sample identifier and must match the IDs used in the trait and covariate files. If you want them to overlap with THGP database measurements, use the `Unique.ID` format `YYYY-MM-DD_BARCODE_ID`.

### 2. Trait file (`traitData`)

Tab-delimited matrix of normalized trait values with **traits in rows** and **samples in columns**:

- The header row is the sample IDs; the first column header is `Unique.ID`.
- Each subsequent row is one trait, identified by its row label (this label becomes the `{trait}` wildcard and the output sub-directory name).
- Values are normalized trait values. Use exactly `NA` or blank (`""` in R) for missing / incorrect / extreme-outlier values.

> **Minimum sample size:** Each trait must have at least **500 non-NA values** among the genotyped samples. Traits below this threshold are skipped automatically and a `results/{trait}/{trait}_SKIPPED.log` file is written explaining why (see `workflow/scripts/orderTrait_byGeno.R`).

### 3. Covariate / metadata file (`metadata`)

Tab-delimited file of covariates to include in the model:

- Column 1: `Unique.ID` — must match IDs in the genotype and trait files.
- Columns 2+: one column per covariate (e.g., `Age`, `Sex`). Column names can be anything; all columns after `Unique.ID` are included as covariates.

Genotype PCs are added to these covariates automatically (the number is set by `GENO_PCnum` in the config).

## Configuration

Edit `config/config.yaml` to set your paths and parameters:

```yaml
##############################
# Variables you need to change
##############################
inputBfile:   ""   # path + PLINK prefix (bfile) to genotype files (no .bed/.bim/.fam)
traitData:    ""   # path to normalized trait matrix (traits in rows, samples in columns)
metadata: ""   # path to covariate/metadata file (Unique.ID + covariate columns)
GENO_PCnum:   5    # number of genotype PCs to include as covariates

##############################
# Variables you *might* need to change
##############################
MAF:  0.05   # minor allele frequency threshold
geno: 0.5    # max per-SNP missingness
mind: 0.95   # minimum per-individual call rate
```

> **Note:** `workflow/Snakefile` reads its config via `configfile: "path/to/config.yaml"`. Update this line to point at your config (e.g. `config/config.yaml`) before running.

## Running the Pipeline

### Dry run (check the workflow)

```bash
snakemake -n --snakefile workflow/Snakefile
```

### Local execution

```bash
snakemake --snakefile workflow/Snakefile --cores 8 --use-conda
```

### SLURM cluster execution

A ready-to-use SLURM profile is provided in `profiles/slurm/`:

```bash
snakemake --snakefile workflow/Snakefile --workflow-profile profiles/slurm
```

The profile sets per-rule memory/runtime resources, enables `--use-conda`, allows up to 1000 jobs, retries failed jobs twice, and groups the lightweight setup rules (`get_PCs`, `get_cov`, `get_plink`) to reduce job overhead.

## Workflow Steps

| Rule | Description |
|------|-------------|
| `get_trait_file` (checkpoint) | Extracts a single trait's values from the trait matrix into `results/{trait}/{trait}_data.tsv` |
| `plink_filter_overlap` | Applies MAF/`geno`/`mind` filters and produces filtered and LD-pruned genotype sets |
| `order_traits` | Orders the trait vector to match the genotype `.fam` ordering and enforces the 500-sample minimum (`orderTrait_byGeno.R`) |
| `get_geno_PCs` | Runs `plink --pca` on the pruned genotypes |
| `get_final_cov` | Combines genotype PCs with user covariates into the GEMMA covariate matrix (`getFinalCovs.R`) |
| `get_LOCO_GRMs` | Builds one standardized kinship matrix per chromosome, each excluding that chromosome |
| `get_SNP_lists` | Splits SNP IDs into per-chromosome lists |
| `get_GEMMA_files` | Assembles per-trait PLINK files with the phenotype written into the `.fam` |
| `run_ME` | Runs `gemma -lmm 4` per trait per chromosome using the matching LOCO kinship matrix, covariates, and SNP list |

## Output Files

### Main results

| File | Description |
|------|-------------|
| `results/{trait}/{trait}_chr{chromosome}_LMM.assoc.txt` | GEMMA LMM association results for one trait on one chromosome |
| `results/{trait}/{trait}_SKIPPED.log` | Written only when a trait has fewer than 500 non-NA samples |

### Intermediate files

| File | Description |
|------|-------------|
| `results/final_geno.*` | Filtered genotype files |
| `results/final_geno_pruned.*` | LD-pruned genotypes for PCA/kinship |
| `results/final_geno_PCA.eigenvec` / `.eigenval` | Genotype principal components |
| `results/nochr{1-22}_GRM.sXX.txt` | LOCO kinship matrices |
| `results/SNPlist_chr{1-22}.txt` | Per-chromosome SNP lists |
| `results/final_GEMMA_cov_W.tsv` | Covariate matrix (PCs + user covariates) |
| `results/{trait}/{trait}_data.tsv` | Extracted per-trait values |
| `results/{trait}/{trait}_genoOverlap.tsv` | Trait values ordered to match genotype samples |

### GEMMA output columns

The `.assoc.txt` files contain standard GEMMA univariate LMM output, including:

`chr`, `rs`, `ps`, `n_miss`, `allele1`, `allele0`, `af`, `beta`, `se`, `logl_H1`, `l_remle`, and the `-lmm 4` test statistics (`p_wald`, `p_lrt`, `p_score`).

## LOCO Strategy

This pipeline uses **Leave-One-Chromosome-Out (LOCO)** kinship matrices to avoid proximal contamination:

- For each chromosome (1–22), a kinship matrix is computed from SNPs on **all other** chromosomes (`plink --not-chr` followed by `gemma -gk 2`).
- When testing SNPs on chromosome *N*, the matrix `nochr{N}_GRM.sXX.txt` is used.
- This prevents the SNPs being tested from also appearing in the kinship (random-effects) term.

Each chromosome's GRM is computed as a separate job for parallelization.

## Troubleshooting

**Sample ID mismatches.** Ensure `Unique.ID` values are consistent across the trait matrix header, the covariate file's first column, and the genotype `.fam` IDs.

**Trait skipped / pipeline stops for a trait.** A trait with fewer than 500 non-NA genotyped samples is intentionally skipped; check the corresponding `_SKIPPED.log`.

**Config not found.** Confirm the `configfile:` path in `workflow/Snakefile` points to your config file.

**Memory issues during kinship computation.** LOCO GRM jobs can be memory-intensive; increase the `get_LOCO_GRMs` allocation in `profiles/slurm/config.yaml` or prune the genotypes more aggressively.

## Citation

If you use this pipeline, please cite GEMMA:

> Zhou X, Stephens M. Genome-wide efficient mixed-model analysis for association studies. *Nature Genetics*. 2012;44(7):821–824.

See the [GEMMA GitHub repository](https://github.com/genetics-statistics/GEMMA) for the latest citation information.
