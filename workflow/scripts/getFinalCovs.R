# get the final covariate file with an intercept, the top genotype PCs, and the input covariates
# rows are in .fam order, as GEMMA requires; missing covariates are written as NA and GEMMA drops those samples

library(tidyverse)

# read in the genotype IDs (.fam order defines the row order of the output)
IDs <- read.table( snakemake@input[["genoIDs"]], sep = " ", header = F)

# read in the genotype PCs (plink --pca eigenvec: FID, IID, PC1, PC2, ...)
PCnum <- snakemake@params[["genoPCnum"]]
genoPCs <- read.table( snakemake@input[["genoPCs"]], sep = " ", header = F)

if ( PCnum > ncol( genoPCs) - 2) {
    stop( paste0( "GENO_PCnum = ", PCnum, " but only ", ncol( genoPCs) - 2, " PCs were computed"))
}

# start with a single intercept column
finalDF <- data.frame( intercept = rep( 1, nrow( IDs)))

# add the top PCs, matched to .fam order
if ( PCnum > 0) {
    pcRows <- match( IDs$V1, genoPCs$V1)
    if ( any( is.na( pcRows))) {
        stop( "some .fam samples are missing from the PCA eigenvec file")
    }
    finalDF <- cbind( finalDF, genoPCs[ pcRows, 3:( PCnum + 2), drop = F])
}

# read in the metadata with the remaining covariates
THGPmeta <- read.table( snakemake@input[["metadat"]], sep = "\t", header = T, check.names = F) %>%
    filter( Unique.ID %in% IDs$V1)

if ( any( duplicated( THGPmeta$Unique.ID))) {
    stop( "duplicated Unique.ID values in the metadata file")
}

# add the additional covariates, matched to .fam order (samples without metadata get NA)
nCov <- ncol( THGPmeta) - 1 # Unique.ID column
if ( nCov > 0) {
    metaRows <- match( IDs$V1, THGPmeta$Unique.ID)
    finalDF <- cbind( finalDF, THGPmeta[ metaRows, -1, drop = F])
}

write.table( finalDF, snakemake@output[["fullCov"]], sep = "\t", quote = F, col.names = F, row.names = F)
