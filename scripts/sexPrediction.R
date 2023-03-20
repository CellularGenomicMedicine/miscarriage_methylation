###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: predict sample sexes based on methylation profiles (for sample QC)

# input: raw beta values and detection p-values (.csv)

# output:  predicted sex per sample

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(sest))
suppressPackageStartupMessages(library(tibble))

# load the data
betas = fread(BetasFilePath)
pvals = fread(PValsFilePath)
ann = fread(annotationFilePath)

# correctly format the betas and pvals for the sest tool
betas = betas %>% column_to_rownames("V1") %>% as.matrix()
pvals = pvals %>% column_to_rownames("V1") %>% as.matrix()

# make the naming column in the annotation file to be able to join the data
ann$Basename = paste0(ann$Sentrix_ID, "_", ann$Sentrix_Position)

# conduct sex estimation
sex = estimateSex(beta.value = betas, detecP = pvals)
ann$Gender = "F" #all samples were recorded as female

# compare the annotated gender to the predicted gender
ann = sex$test %>% rownames_to_column("Basename") %>% 
  select(Basename, predicted.X, predicted.Y, predicted, X.PC1, Y.PC1) %>%
  full_join(ann, ., by = "Basename")
  
write.csv(ann, annotationOutFilePath.csv, row.names = F)
