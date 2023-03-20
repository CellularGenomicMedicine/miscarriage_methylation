###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: combat batch correction to correct the effect associated with data from 2 different runs

# input: preprocessed beta values (.csv) & relevant annotation files (.csv)

# output: combat corrected beta values per site per sample (.csv)

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(tibble))

# load the annotation file
ann = fread(AnnotationFilePath.csv)

# join the betas of the 2 separate runs (load all the data and filter both to contain the same sites, then join)

betas1 = fread(Run1BetasFilePath) 

betas2 = fread(Run2BetasFilePath) %>%
        filter(V1 %in% betas1$V1)
betas1 = betas1 %>% filter(V1 %in% betas2$V1)

betas = full_join(betas1, betas2, by = "V1") %>%
        column_to_rownames("V1")

# format the data for combat
betas = as.matrix(betas)

# apply combat
betas = ComBat(betas, batch = ann$run)

# save the output
write.csv(betas, file = file.path(report.dir, "betasCombat.csv"), row.names = T)

