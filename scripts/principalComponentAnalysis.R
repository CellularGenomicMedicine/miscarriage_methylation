###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: conduct PCA and associate PCs with sample characteristics

# input: fully pre-processed beta values (high quality, normalised and batch corrected), sample annotation file

# output:  PCA plot, associataions heatmap

###########################################################################################################################

#load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(jmuOutlier))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(tidyr))

# load the sample annotation
ann = fread(annotationFile)

# load the processed beta values
betas = fread(processedBetasFile)

# prepare the data for PCA (omit any sites containing missing values)
betas.pca <- betas %>% na.omit() %>%
  as.matrix() %>%
  t()
  
 # conduct a principal component analysis
pca <- prcomp(betas.pca, center = T, scale. = F)

# to see % of variance explained by each PC
summary(pca)

# join the PC coordinates to the sample annotation
coords <- as.data.frame(pca$x) %>% rownames_to_column("Sample_Name") %>%
  full_join(ann, .[ , 1:15], by = "Sample_Name")
  
# plot the results

pdf(file = "outPutFile.pdf", width = 10, height = 6)

plot = coords %>% ggplot(aes(x = PC1, y = PC2, colour = individual, shape = TissueType)) + 
  geom_point(size = 2) +
  scale_shape_manual(values = c("CV" = 15, "EM" = 1), name = "Tissue type") +
  #scale_colour_manual(values = c("CV" = "#799141", "EM" = "#E67A74"), name = "Tissue type") + 
  stat_ellipse(aes(group = TissueType), level = 0.90) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  xlab("PC1 (50.1%)") +
  ylab("PC2 (15.2%)") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(family = "Helvetica"))
        
dev.off()
  
#### Calculate associations between the PCs and available sample characteristics

# run (binary categorical = Wilcoxon rank)
run = list()

for(n in 1:5){
  pc = paste0("PC", n)
  x = coords %>% filter(run == 1) %>% select(all_of(pc)) 
  x = x %>% as.matrix() %>% .[,1]  
  y = coords %>% filter(run == 2) %>% select(all_of(pc))
  y = y %>% as.matrix() %>% .[,1] 
  test = wilcox.test(x = x, y = y, alternative = "two.sided")  
  run[[pc]] = test$p.value  
}

# gestational age
ga = list()
cor.ga = list()

for(n in 1:5){
  pc = paste0("PC", n) 
  data = coords %>% select(Gest_age_week, all_of(pc)) %>% na.omit()
  test = perm.cor.test(x = unlist(data[,1]), y = unlist(data[,2]), alternative = "two.sided", 
                       method = "pearson", num.sim = 10000)
  ga[[pc]] = test$p.value
  cor.ga[[pc]] = cor(data[,1], data[,2], method = "pearson")
}

# tissue type
tissue = list()

for(n in 1:5){  
  pc = paste0("PC", n)
  x = coords %>% filter(TissueType == "CV") %>% select(all_of(pc)) 
  x = x %>% as.matrix() %>% .[,1] 
  y = coords %>% filter(TissueType == "EM") %>% select(all_of(pc)) 
  y = y %>% as.matrix() %>% .[,1]
  test = wilcox.test(x = x, y = y, alternative = "two.sided")
  tissue[[pc]] = test$p.value
}


#Cell composition

Stromal = list()
cor.Stromal = list()

for(n in 1:5){
  pc = paste0("PC", n)
  data = coords %>% select(Stromal, all_of(pc))
  test = perm.cor.test(x = unlist(data[,1]), y = unlist(data[,2]), alternative = "two.sided", method = "pearson", num.sim = 10000)
  Stromal[[pc]] = test$p.value
  cor.Stromal[[pc]] = cor(data[,1], data[,2], method = "pearson")
}

Hofbauer = list()
cor.Hofbauer = list()

for(n in 1:5){
  pc = paste0("PC", n)
  data = coords %>% select(Hofbauer, all_of(pc))
  test = perm.cor.test(x = unlist(data[,1]), y = unlist(data[,2]), alternative = "two.sided", method = "pearson", num.sim = 10000)
  Hofbauer[[pc]] = test$p.value
  cor.Hofbauer[[pc]] = cor(data[,1], data[,2], method = "pearson")
}

Endothelial = list()
cor.Endothelial = list()

for(n in 1:5){
  pc = paste0("PC", n)
  data = coords %>% select(Endothelial, all_of(pc))
  test = perm.cor.test(x = unlist(data[,1]), y = unlist(data[,2]), alternative = "two.sided", method = "pearson", num.sim = 10000)
  Endothelial[[pc]] = test$p.value
  cor.Endothelial[[pc]] = cor(data[,1], data[,2], method = "pearson")
}

nRBC = list()
cor.nRBC = list()

for(n in 1:5){
  pc = paste0("PC", n)
  data = coords %>% select(nRBC, all_of(pc)) 
  test = perm.cor.test(x = unlist(data[,1]), y = unlist(data[,2]), alternative = "two.sided", method = "pearson", num.sim = 10000)
  nRBC[[pc]] = test$p.value
  cor.nRBC[[pc]] = cor(data[,1], data[,2], method = "pearson")
}

Syncytiotrophoblast = list()
cor.Syncytiotrophoblast = list()

for(n in 1:5){
  pc = paste0("PC", n)
  data = coords %>% select(Syncytiotrophoblast, all_of(pc))
  test = perm.cor.test(x = unlist(data[,1]), y = unlist(data[,2]), alternative = "two.sided", method = "pearson", num.sim = 10000)
  Syncytiotrophoblast[[pc]] = test$p.value
  cor.Syncytiotrophoblast[[pc]] = cor(data[,1], data[,2], method = "pearson")
}

## Combine the output from the statistical testing
associations = data.frame(row.names = c("PC1", "PC2", "PC3", "PC4", "PC5")) %>%
  rownames_to_column("Principal_component") %>%
  mutate(run = unlist(run),
         GestationalAge = unlist(ga),
         TissueType = unlist(tissue),
         Stromal = unlist(Stromal),
         Hofbauer = unlist(Hofbauer),
         Endothelial = unlist(Endothelial),
         nRBC = unlist(nRBC),
         Syncytiotrophoblast = unlist(Syncytiotrophoblast))
         
## Re-format the data for plotting
to.plot = associations %>% tibble::column_to_rownames("Principal_component") %>% t()

t = to.plot %>% as.data.frame() %>% rownames_to_column("variable")

t2 = gather(t, key = "PC", value = "pvals", PC1:PC5)
t2$variable = factor(t2$variable, levels = rev(rownames(to.plot)), ordered = T)

t2 = t2 %>% mutate(logpvals = abs(log10(pvals)), 
                   text = pvals <= 0.05, 
                   pval.text = round(pvals, digits = 3))
                   
## Generate the heatmap
pdf(file = "outPutFile.pdf", width = 6, height = 3)
t2 %>% ggplot(aes(x = PC, y = variable, fill = logpvals)) + 
  geom_tile(colour = "white") + 
  scale_fill_gradientn(colours = c("white", "lightblue", "dodgerblue2"), 
                       values = scales::rescale(c(0, -log10(0.05), 8)), limits = c(0,5), oob = scales::squish) +
  geom_text(data = function(x){filter(x, text)}, aes(label = pval.text), size = 2, angle = 20) +
  theme_classic() +
  labs(fill = "-log10 \n p-value") +
  scale_y_discrete(labels = c("GestationalAge" = "Gestational age", "TissueType" = "Tissue type", "Stromal" = "Stromal cells",
                              "Hofbauer" = "Hofbauer cells", "Endothelial" = "Endothelial cells", "nRBC" = "nRBC",
                              "Syncytiotrophoblast" = "Syncytiotrophoblast cells")) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        text = element_text(family = "Helvetica", size = 7),
        plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        axis.text.y = element_text(size = 7)
        ) +
  theme(text = element_text(family = "Helvetica"))
dev.off()
