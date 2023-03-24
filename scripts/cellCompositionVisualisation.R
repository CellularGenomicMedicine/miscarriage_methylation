###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualise the predicted cell compositions of the samples

# input: sample annotation information with the appended cell composition estimates

# output: box plots per cell type arranged in a grid and saved as a .pdf

###########################################################################################################################


# load packages

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

# load the annotation data (with the appended predicted cell compositions)
ann = read.csv

# set the desired colour scheme
colours = data.frame(cellType = c("Syncytiotrophoblast", "Stromal", "nRBC", "Hofbauer", "Endothelial"),
                     shading = c("#f3c300", "#875692", "#f38400", "#a1caf1", "#be0032"))
                     
# generate a plot per cell type (Trophoblasts are not included as they were not predicted to be present in our samples)

plots = list()

for(c in unique(ann %>% filter(cellType != "Trophoblasts") %>% .$cellType)){
  
  cellColour = colours %>% filter(cellType == c) %>% .$shading
  
  plot = ann %>% filter(cellType == c) %>%
    ggplot(aes(x = TissueType, y = proportion)) +
    #geom_flat_violin(aes(fill = TissueType, colour = TissueType), 
    #                 position = position_nudge(x = 0, y = 0), adjust = 0.75, trim = FALSE) +
    geom_boxplot(aes(group = TissueType), width = 0.8, outlier.shape = NA, colour = "black", lwd = 0.4, fill = cellColour) +
    geom_point(aes(group = TissueType), size = 1) +
    theme_classic() +
    ggtitle(c) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(family = "Helvetica", size = 11))
  
  plots[[c]] = plot
  
}

# make a .pdf file containing a grid of the plots

pdf(file = "outPutFilePath.pdf", width = 12, height = 3)
plot_grid(plots[[1]] + theme(axis.title.x = element_blank()), 
          plots[[2]] + theme(axis.title.y = element_blank(),
                             axis.title.x = element_blank()), 
          
          plots[[3]] + theme(axis.title.y = element_blank(),
                             axis.title.x = element_blank()),
          
          plots[[4]] + theme(axis.title.y = element_blank(),
                           axis.title.x = element_blank()), 
          
          plots[[5]] + theme(axis.title.y = element_blank(),
                             axis.title.x = element_blank()),
          nrow = 1, rel_widths = c(1.1, 1, 1),
          labels = c("a)", "b)", "c)", "d)", "e)"),
          label_size = 11, label_fontfamily = "Helvetica")
dev.off()
