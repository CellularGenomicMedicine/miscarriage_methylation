###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: cell composition prediction based on reference data

# input: preprocessed (methylumiNoob) beta values (.csv), sample annotation (.csv), reference data (see supporting data
# directory)

# Supporting data reference: Yuan V, Hui D, Yin Y, Penaherrera MS, Beristain AG, Robinson WP. Cell-specific characterisation
# of the placental methylome. BMC Genomics. 2021.

# output:  cell composition prediction (cytotrophoblasts, syncytiotrophoblasts, nucleated red blood cells, Hofbauer cells,
# endothelial cells, stromal cells)

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))

# load the data
annotation = fread(AnnotationFilePath.csv)
betas = fread(BetasFilePath)
ref = fread(referenceData) %>%
        column_to_rownames("cpg") %>%
        as.matrix()
        
# filter the beta values to only contain the sites in the reference data
betas.ref = betas %>% filter(rownames(.) %in% rownames(ref)) %>%
        as.matrix()
        
# calculate the cell estimates using the Houseman method from minfi
cells = minfi:::projectCellType(
        betas.ref,
        ref)
        
# add the cell compositions to the sample annotation
cells = cells %>% as.data.frame() %>%  rownames_to_column("Basename")
full = full_join(annotation, cells, by = "Basename")    

# save the updated sample sheet

write.csv(full, outputFileLocation.csv, row.names = F)
