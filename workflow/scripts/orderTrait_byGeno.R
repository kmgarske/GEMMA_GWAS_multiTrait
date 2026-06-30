# get the RNA expression data in the same order as the genotyping file

library(tidyverse)

# read in the genotyped samples
geno <- read.table( snakemake@input[["genoIDs"]], sep = " ", header = F)

# read in the gene expression data
trait <- read.table( snakemake@input[["traitData"]], sep = "\t", header = T, check.names = F) %>% filter( Unique.ID %in% geno$V1)

finalDF <- data.frame( Unique.ID = geno$V1) %>% left_join( trait, by = "Unique.ID") %>% arrange( match( Unique.ID, geno$V1))

# check for minimum sample size (non-NA values)
traitCol <- finalDF[[2]]
nNonNA <- sum( !is.na( traitCol))
minSamples <- 500

if ( nNonNA < minSamples) {
    # write a skip log file so the user knows this was intentional
    skipLogFile <- gsub( "_genoOverlap.tsv$", "_SKIPPED.log", snakemake@output[["finalDat"]])
    skipMsg <- paste0( "TRAIT SKIPPED: Insufficient samples\n",
                       "Non-NA samples found: ", nNonNA, "\n",
                       "Minimum required: ", minSamples, "\n",
                       "Date: ", Sys.time(), "\n")
    writeLines( skipMsg, skipLogFile)
    message( skipMsg)
    stop( paste0( "Insufficient samples: only ", nNonNA, " non-NA values found (minimum required: ", minSamples, "). ",
                  "Stopping pipeline for this trait. See: ", skipLogFile))
}

# save everything
write.table( finalDF %>% select( 2), snakemake@output[["finalDat"]], sep = "\t", quote = F, col.names = F, row.names = F)
