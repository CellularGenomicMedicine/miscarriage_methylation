###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualise the methylation value distribution across samples

# input: processed beta values, sample annotation file

# output: density plot of beta values per sample / condition

###########################################################################################################################

#load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tibble))

# set data directories
annotation.file = "sampleAnnotationFile.csv"

betas.file = "processedBetasFile.csv"

# load the data files

annotation <- fread(annotation.file) 

betas <- fread(betas.file) %>% 
  select(ID, all_of(annotation$Sample_Name)) %>%
  column_to_rownames("ID")

# plot the beta vlue distribution of each sample, coloured by tissue type
betas %>% ggplot() +
  geom_density(aes(x = PL2712_EM3700), colour = "#E67A74", fill = NA) +
  geom_density(aes(x = PL2712_CV3700), colour = "#799141", fill = NA) +
  geom_density(aes(x = PL2617_CV3627), colour = "#799141", fill = NA) +
  geom_density(aes(x = PL2717_EM3704), colour = "#E67A74", fill = NA) +
  geom_density(aes(x = PL2717_CV3704), colour = "#799141", fill = NA) +
  geom_density(aes(x = PL2682_EM3669), colour = "#E67A74", fill = NA) +
  geom_density(aes(x = PL2682_CV3669), colour = "#799141", fill = NA) +
  geom_density(aes(x = PL245_EM1081), colour = "#E67A74", fill = NA) +
  geom_density(aes(x = PL245_CV1081), colour = "#799141", fill = NA) +
  geom_density(aes(x = PL2019_EM3509), colour = "#E67A74", fill = NA) +
  geom_density(aes(x = PL2137_CV3531), colour = "#799141", fill = NA) +
  geom_density(aes(x = PL2409_EM3584), colour = "#E67A74", fill = NA) +
  geom_density(aes(x = PL2409_CV3584), colour = "#799141", fill = NA) +
  theme_classic()

# re-format the data to facilitate average plotting per condition
betas2 = betas %>%
   t()

betas3 = betas2 %>% as.data.frame() %>% rownames_to_column("Sample_Name")

ann = annotation %>% select(Sample_Name, TissueType) 

betas3 = betas3 %>% mutate(Group = ann[Sample_Name == betas3$Sample_Name, TissueType])

long = tidyr::gather(betas3, ID, beta, 2:687216)

# plot the total methylation value distribution per tissue type
long %>% 
  ggplot(aes(x = beta, colour = Group)) + 
  geom_density(fill = NA, alpha = 0.7) +
  scale_colour_manual(values = c("CV" = "#799141", "EM" = "#E67A74"), name = "Tissue Type") +
  theme_classic()


###################### Calculate & visualise the mean methylation per sample #########################

## calculate the average methylation per sample

means = rowMeans(betas2, na.rm = T)

annotation$meanMeth = means

annotation %>% ggplot(aes(x = TissueType, y = meanMeth)) +
  geom_violin(aes(fill = TissueType)) +
  scale_fill_manual(values = c("CV" = "#799141", "EM" = "#E67A74")) +
  geom_boxplot(width = 0.1) +
  geom_point(colour = "grey") +
  theme_classic()
