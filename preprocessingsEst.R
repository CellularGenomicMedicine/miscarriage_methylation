###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: preprocessing of raw .idat files to carry out sex prediction with SEst

# input: .idat files and associated sample sheets (.csv)

# output:  raw methylation beta value per site per sample, detection p-values per site per sample (format: .csv)

###########################################################################################################################

# load the required packages
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

# set the data directories
data.dir <- DataFilePath
sample.annotation <- AnnotationFilePath.csv
report.dir <- OutputFilePath

# read the raw data files using minfi

targets <- fread(sample.annotation) %>%
  mutate(Basename = paste(Sentrix_ID, Sentrix_Position , sep = "_"))

RGSet <- read.metharray.exp(base = data.dir, targets = targets)

# proces the RGSet to an MSet  which has beta values

MSet <- preprocessRaw(RGSet)

# calculate beta values by converting this to a ratioSet object

RSet <- ratioConvert(MSet, what = "beta")

# extract beta values and detection p values for all sites

detPval <- detectionP(RGSet)

betas <- getBeta(RSet)

# write the resultant files as csv files for further processing

write.csv(detPval, file = file.path(report.dir,  "detPval.csv"), row.names = T)

write.csv(betas, file = file.path(report.dir, "betasRaw.csv"), row.names = T)
