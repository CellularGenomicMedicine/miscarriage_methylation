###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: look at methylation at repetitive elements alu and line1 (can be compared to reference data from Yuan et al)

# input: Alu and Line1 sites, processed beta values, sample annotation

# output: box plots of average methylation per repetitive element per tissue type

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

##################### for the alu sites #####################################
# load the data
alu = data.table::fread("aluSites.csv")

betas = data.table::fread("betasFile.csv") %>%
  select(-Chromosome, -Start, -End, -Strand)

ann = data.table::fread("annotationFile.csv") 

# select the relevant sites

betas_alu = betas %>% filter(ID %in% alu$ID,
                             ID %in% refbetas$ID_REF)

# arrange the data for plotting

betas_alu = betas_alu %>% tibble::column_to_rownames("ID") %>%
  t()

summary_alu = data.frame(Sample_Name = rownames(betas_alu),
                            means = rowMeans(betas_alu, na.rm = T))

all = full_join(ann, summary_alu, by = "Sample_Name")                           

# make the summaries suitable for joining

summary_alu = summary_alu %>% select(Sample_Name, TissueType, means)

sum_alu_ref = sum_alu_ref %>% filter(Trimester == "First") %>%
  select(Sample_ID, CellType, means)

colnames(sum_alu_ref) = c("Sample_Name", "TissueType", "means")

all = rbind(summary_alu, sum_alu_ref)

# plot the data
all %>% ggplot() +
  geom_violin(aes(y = means, x = TissueType, colour = TissueType, fill = TissueType)) +
  geom_boxplot(aes(x = TissueType, y = means), width = 0.1) +
  geom_point(aes(x = TissueType, y = means), size = 0.5, colour = "grey") +
  scale_fill_manual(values = c("CV" = "#799141", "EM" = "#E67A74")) +
  scale_colour_manual(values = c("CV" = "#799141", "EM" = "#E67A74")) +
  theme_classic() +
  theme(legend.position = "none")

#########################################################################

# Line1 repetitive elements

# load the data
line = data.table::fread("lineSites.csv")

# select the relevant sites

betas_line = betas %>% filter(ID %in% line$ID,
                             ID %in% refbetas$ID_REF)

# arrange the data for plotting

betas_line = betas_line %>% tibble::column_to_rownames("ID") %>%
  t()

summary_line = data.frame(Sample_Name = rownames(betas_line),
                         means = rowMeans(betas_line, na.rm = T))

all = full_join(ann, summary_line, by = "Sample_Name")                           

# plot the data
all_line %>% ggplot() +
  geom_violin(aes(y = means, x = TissueType, colour = TissueType, fill = TissueType)) +
  geom_boxplot(aes(x = TissueType, y = means), width = 0.1) +
  geom_point(aes(x = TissueType, y = means), size = 0.5, colour = "grey") +
  scale_fill_manual(values = c("CV" = "#799141", "EM" = "#E67A74")) +
  scale_colour_manual(values = c("CV" = "#799141", "EM" = "#E67A74")) +
  theme_classic() +
  theme(legend.position = "none")

