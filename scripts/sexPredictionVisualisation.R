###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualisation of sex prediction output

# input: sex prediction from SEst tool

# output:  plot (saved as .RDS file for later formatting / multi-panel figure generation)

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

# load the data
ann = fread(AnnotationFileLocation.csv)

# plot the data
plot = ann %>% ggplot(aes(x = X.PC1, y = Y.PC1, colour = predicted, shape = Gender)) +
  geom_point() +
  theme_classic() +
  scale_shape_manual(values = c("F" = 16, "M" = 17), labels = c("Female", "Male")) +
  scale_color_manual(values = c("F" = "plum2", "M" = "skyblue2", "N" = "grey"), labels = c("Female", "Male", "Not specified")) +
  labs(shape = "Recorded sex", colour = "Predicted Sex") +
  scale_x_continuous(limits = c(-0.3, 0.3), breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.5, 0.7), breaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6), expand = c(0, 0.05))

saveRDS(plot, file = outFilePath/sEST.rds")
