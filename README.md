# GEMMA_GWAS
Run a GEMMA GWAS across all traits in the Turkana Health and Genomics Project (THGP) in an input trait matrix with input plink files

Implements leave-one-chromosome-out (LOCO) method.

## Input file formats

### genotype files
PLINK format bed, bim, fam for all samples that have been genotyped, with names as Unique.ID in YYYY-MM-DD_BARCODE_ID format if you want them to overlap with the THGP database measurements

### trait file
tab-delimited
column names: Unique.ID in YYYY-MM-DD_BARCODE_ID format from the THGP database

row names: trait

values: normalized trait values, exactly NA or blank ("" in R) if missing/incorrect/extreme outlier/etc

### covariate file
tab-delimited

column 1 name: Unique.ID

column 1 values: Unique.ID in YYYY-MM-DD_BARCODE_ID format from the THGP database

column 2 name: <covariate1Name> e.g., "Age" - can be anything

column 2 values: value for <covariate1> 

column 3 name: <covariate2Name> e.g., "Sex" - can be anything

column 3 values: value for <covariate2>

...and so on for all covariates you want to include in your model. genotype PCs will also be included (see config/config.yaml file)
