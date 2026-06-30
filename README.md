# GEMMA_GWAS
Run a GEMMA GWAS across all traits in an input trait matrix with input plink files and metadata.

Implements leave-one-chromosome-out (LOCO) method.

## Input file formats

### genotype files
PLINK format bed, bim, fam for all samples that have been genotyped, with names as Unique.ID in YYYY-MM-DD_BARCODE_ID format if you want them to overlap with the THGP database measurements

### trait file
tab-delimited
column names: sampleID should match IDs in plink and metadata

row names: trait

values: normalized trait values, exactly NA or blank ("" in R) if missing/incorrect/extreme outlier/etc

### covariate file
tab-delimited

column 1 name: sampleID

column 1 values: sampleID should match IDs in plink and metadata

column 2 name: <covariate1Name> e.g., "Age" - can be anything

column 2 values: value for <covariate1> 

column 3 name: <covariate2Name> e.g., "Sex" - can be anything

column 3 values: value for <covariate2>

...and so on for all covariates you want to include in your model. genotype PCs will also be included (see config/config.yaml file)

## Citation

If you use this pipeline, please cite:
- Zhou X, Stephens M. Genome-wide efficient mixed-model analysis for association studies. *Nature Genetics*. 2012;44(7):821-824.
