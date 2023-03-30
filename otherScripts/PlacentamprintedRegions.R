###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualise the methylation values of sites specifically methylated in the placenta or all tissues

# input: site information (pladental and generally imprinted sites), processed beta values, sample annotation file

# output: density plots showing the distribution of beta values at imprinted sites, per tisuse or per sample

###########################################################################################################################

# load the required packages

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# load the data

sites = data.table::fread("placentaImprintedSites.csv")

placenta_sites = sites %>% filter(tissue_specificity == "placental-specific") %>% .$ID
other_sites = sites %>% filter(tissue_specificity == "other") %>% .$ID

betas = data.table::fread("processedBetasFile.csv") %>%
  select(-Chromosome, -Start, -End, -Strand)

ann = data.table::fread("sampleAnnotationFile.csv") %>%
  select(-V1) 

# select the relevant sites
betas_sites = betas %>% filter(ID %in% sites$ID)

# merge the sites with the information about the imprinting
betas_sites = full_join(betas_sites, sites %>% filter(ID %in% betas_sites$ID), by = "ID")

# remove duplicated rows

betas_sites = distinct(betas_sites)

# calculate the methylation level of the relevant sites
## imprinted is 0.25 - 0.75

imprinted = betas_sites

imprinted$PL2712_EM3700 = if_else(imprinted$PL2712_EM3700 > 0.25 & imprinted$PL2019_EM3509 < 0.75, "IMP", "no")
imprinted$PL2712_CV3700 = if_else(imprinted$PL2712_CV3700 > 0.25 & imprinted$PL2712_CV3700 < 0.75, "IMP", "no")
imprinted$PL2617_CV3627 = if_else(imprinted$PL2617_CV3627 > 0.25 & imprinted$PL2617_CV3627 < 0.75, "IMP", "no")
imprinted$PL2717_EM3704 = if_else(imprinted$PL2717_EM3704 > 0.25 & imprinted$PL2717_EM3704 < 0.75, "IMP", "no")
imprinted$PL2717_CV3704 = if_else(imprinted$PL2717_CV3704 > 0.25 & imprinted$PL2717_CV3704 < 0.75, "IMP", "no")
imprinted$PL2682_EM3669 = if_else(imprinted$PL2682_EM3669 > 0.25 & imprinted$PL2682_EM3669 < 0.75, "IMP", "no")
imprinted$PL2682_CV3669 = if_else(imprinted$PL2682_CV3669 > 0.25 & imprinted$PL2682_CV3669 < 0.75, "IMP", "no")
imprinted$PL245_EM1081 = if_else(imprinted$PL245_EM1081 > 0.25 & imprinted$PL245_EM1081 < 0.75, "IMP", "no")
imprinted$PL245_CV1081 = if_else(imprinted$PL245_CV1081 > 0.25 & imprinted$PL245_CV1081 < 0.75, "IMP", "no")
imprinted$PL2019_EM3509 = if_else(imprinted$PL2019_EM3509 > 0.25 & imprinted$PL2019_EM3509 < 0.75, "IMP", "no")
imprinted$PL2137_CV3531 = if_else(imprinted$PL2137_CV3531 > 0.25 & imprinted$PL2137_CV3531 < 0.75, "IMP", "no")
imprinted$PL2409_EM3584 = if_else(imprinted$PL2409_EM3584 > 0.25 & imprinted$PL2409_EM3584 < 0.75, "IMP", "no")
imprinted$PL2409_CV3584 = if_else(imprinted$PL2409_CV3584 > 0.25 & imprinted$PL2409_CV3584 < 0.75, "IMP", "no")

# look at placenta-specific imprints
placenta = imprinted %>% filter(tissue_specificity == "placental-specific")

# look at the generally imprinted methylation sites
other = imprinted %>% filter(tissue_specificity != "placental-specific") 

# reformat the data for plotting
betas2 = betas_sites %>% select(ID, all_of(ann$Sample_Name)) %>%
  tibble::column_to_rownames("ID") %>% 
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("Sample_Name")

betas2 = ann %>% select(Sample_Name, TissueType) %>%
  full_join(., betas2, by = "Sample_Name")

long = tidyr::gather(betas2, ID, beta, 3:ncol(betas2))

# plot the placenta-specific values

long_placenta = long %>% filter(ID %in% placenta_sites)

ggplot() +
  geom_density(data = long_placenta,
               aes(x = beta, colour = TissueType), size = 1.5) +
  facet_wrap(~TissueType) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "50%", "", "100%"),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7),
                     breaks = c(0, 2, 4, 6),
                     labels = c(0, 2, 4, 6),
                     expand = c(0,0)) +
  xlab("DNAm") +
  ggtitle("Imprinted in placenta") +
  scale_colour_manual(values = c("CV" = "#799141", "EM" = "#E67A74"), name = "Tissue type") 

# plot the sites that are imprinted in many tissues

long_other = long %>% filter(ID %in% other_sites)

ggplot() +
  geom_density(data = long_other,
               aes(x = beta, colour = TissueType), size = 1.5) +
  facet_wrap(~TissueType) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "50%", "", "100%"),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7),
                     breaks = c(0, 2, 4, 6),
                     labels = c(0, 2, 4, 6),
                     expand = c(0,0)) +
  xlab("DNAm") +
  ggtitle("Imprinted in multiple tissues") +
  scale_colour_manual(values = c("CV" = "#799141", "EM" = "#E67A74"), name = "Tissue type") 


## plot the data showing the individual samples (check for sample variability)

# placenta sites
ggplot() +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2712_EM3700"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2712_CV3700"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2617_CV3627"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2717_EM3704"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2717_CV3704"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2682_EM3669"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2682_CV3669"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL245_EM1081"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL245_CV1081"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2019_EM3509"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2137_CV3531"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2409_EM3584"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_placenta %>% filter(Sample_Name == "PL2409_CV3584"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
    facet_wrap(~TissueType) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(0, 1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("", "", "50%", "", "100%"),
                       expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 7),
                       breaks = c(0, 2, 4, 6),
                       labels = c(0, 2, 4, 6),
                       expand = c(0,0)) +
    xlab("DNAm") +
    ggtitle("Imprinted in placenta") 

# other sites
ggplot() +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2712_EM3700"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2712_CV3700"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2617_CV3627"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2717_EM3704"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2717_CV3704"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2682_EM3669"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2682_CV3669"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL245_EM1081"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL245_CV1081"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2019_EM3509"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2137_CV3531"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2409_EM3584"), aes(x = beta), colour = "#E67A74", fill = NA, size = 1) +
  geom_density(data = long_other %>% filter(Sample_Name == "PL2409_CV3584"), aes(x = beta), colour = "#799141", fill = NA, size = 1) +
  facet_wrap(~TissueType) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 1), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("", "", "50%", "", "100%"),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7),
                     breaks = c(0, 2, 4, 6),
                     labels = c(0, 2, 4, 6),
                     expand = c(0,0)) +
  xlab("DNAm") +
  ggtitle("Imprinted in multiple tissues") 
