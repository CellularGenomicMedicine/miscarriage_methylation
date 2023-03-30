###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: look at placenta-specific partially methylated domains

# input: sites relating to placenta-specific partially methylated domains, annotation, processed beta values

# output: plot representing placenta-specific partial methylation

############################################################################################################################ script to visualise placental DMPs

# load packages

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))


# load the data

regions = data.table::fread("placentaPMDsites.csv")

betas = data.table::fread("processedBetasFile.csv")

ann = data.table::fread("SampleAnnotationFile.csv") 

# select the relevant CpG sites

betas = betas %>% filter(ID %in% regions$ID) 

betas2 = betas %>% select(-Chromosome, -Start, -End, -Strand) %>%
  tibble::column_to_rownames("ID") %>%
  t() 

# calculate what percentage of sites fall into each methylation range (<0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8, >0.8)
summary = data.frame(Sample_Name = rownames(betas2),
                     lowest = rowSums(betas2 < 0.2, na.rm = T) / ncol(betas2) * 100,
                     low = rowSums(betas2 >= 0.2 & betas2 < 0.4, na.rm = T) / ncol(betas2) * 100,
                     med = rowSums(betas2 >= 0.4 & betas2 < 0.6, na.rm = T) / ncol(betas2) * 100,
                     high = rowSums(betas2 >= 0.6 & betas2 < 0.8, na.rm = T) / ncol(betas2) * 100,
                     highest = rowSums(betas2 >= 0.8, na.rm = T) / ncol(betas2) * 100)

# attach the sample annotation

summary = full_join(ann, summary, by = "Sample_Name")

# rearrange the table for plotting

long = tidyr::gather(summary, level, percentage, lowest:highest)

long$value = long$level
long$value = gsub("lowest", 10, long$value)
long$value = gsub("low", 30, long$value)
long$value = gsub("med", 50, long$value)
long$value = gsub("highest", 90, long$value)
long$value = gsub("high", 70, long$value)

long$value = as.numeric(long$value)

avg = long %>% group_by(TissueType, value) %>% summarise(mean = mean(percentage))

# plot the results

ggplot() +
  geom_point(data = long, aes(x = value, y = percentage, colour = TissueType), size = 1) +
  geom_line(data = avg, aes(x = value, y = mean, colour = TissueType), size = 1.5) +
  geom_point(data = avg, aes(x = value, y = mean, colour = TissueType), size = 3) +
  scale_y_continuous(limits = c(0, 50), breaks = c(0, 25, 50), labels = c("0%", "25%", "50%"), expand = c(0,0)) +
  scale_x_continuous(limits = c(5,95), breaks = c(10, 30, 50, 70, 90), labels = c("<20", "20-40", "40-60", "60-80", ">80"),
                     expand = c(0,0)) +
  theme_classic() +
  xlab("DNAm interval (%)") +
  ylab("Percentage of CpGs") +
  scale_colour_manual(values = c("CV" = "#799141", "EM" = "#E67A74"), name = "Tissue type")
