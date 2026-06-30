# get the final covariate file with input covariates and the top genotype PCs

library(tidyverse)

# read in the genotype PCs
PCnum <- snakemake@params[["genoPCnum"]]
finalCol <- PCnum + 2
genoPCs <- read.table( snakemake@input[["genoPCs"]], sep = " ", header = F) %>% mutate( addCol = 1) 

# read in the genotype IDs
IDs <- read.table( snakemake@input[["genoIDs"]], sep = " ", header = F)

# read in the THGP metadata with the remaining covariates
# THGPmeta <- read.table( snakemake@input[["THGPcov"]], sep = "\t", header = T) %>% filter( Unique.ID %in% genoPCs$V1) %>% arrange( match( Unique.ID, genoPCs$V1))
THGPmeta <- read.table( snakemake@input[["THGPcov"]], sep = "\t", header = T) %>% filter( Unique.ID %in% genoPCs$V1) 

genoPCs <- genoPCs[ , c( ncol( genoPCs), 3:finalCol)]

# add the additional covariates to the PCs
nCov <- ncol( THGPmeta) - 1 # Unique.ID column
if ( nCov > 0) {

    # genoPCs <- cbind( genoPCs, THGPmeta[ , seq( 2, ( 1 + nCov))])
    finalDF <- left_join( IDs %>% rename( Unique.ID = V1), THGPmeta, by = "Unique.ID") %>% arrange( match( Unique.ID, genoPCs$V1))
    finalDF <- cbind( genoPCs, finalDF[ , seq( 7, ( 6 + nCov))])
}

# write.table( genoPCs, snakemake@output[["fullCov"]], sep = "\t", quote = F, col.names = F, row.names = F)
write.table( finalDF, snakemake@output[["fullCov"]], sep = "\t", quote = F, col.names = F, row.names = F)