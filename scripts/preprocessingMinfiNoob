###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: preprocessing of raw .idat files using methylumiNoob (recommended for cell type deconvolution)

# input: .idat files from each run (1 or 2) processed separately, associated sample sheets (.csv)

# output:  methylation beta value per site per sample (format: .csv)

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

#set file directories
idat.dir <- DataFilePath
sample.annotation <- AnnotationFilePath.csv
report.dir <- OutputFilePath

# load the annotation file

annotation <- read.csv(annotation.file)

annotation$Basename = paste0(annotation$Sentrix_ID, "_", annotation$Sentrix_Position)

# load the idat files

RGSet = read.metharray.exp(data.dir, targets = annotation)

# process the data using preprocessNoob

message("preprocessing")

MSet = preprocessNoob(RGSet)

message("converting to ratio set")

ratioSet = ratioConvert(MSet)

# extract only the beta values

betas = getBeta(ratioSet)

# write the betas file out

write.csv(betas, file = file.path(report.dir, "betas.csv"), row.names = T)
