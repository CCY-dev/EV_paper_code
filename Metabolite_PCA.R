rm(list = ls())
setwd("EV_CKD/20220207_version/Analysis/Other_analysis/")
library(tidyverse)
library(phytools)
library(vegan)
library(ade4)
library(readxl)
library(ape)
library(ggpubr)
library(factoextra)

raw <- read_xlsx("../../Data/Database_CKD-EVs_polished_noRed_noBlue.xlsx", 
                 col_names = T)
metabo <- raw[ , c(1, 2, 61:91)]

# Replace "NA"s with NAs and remove the row with NAs
metabo[metabo == "NA"] <- NA
sum(is.na(metabo))
metabo <- metabo %>% 
  drop_na()

### Do PCA
# Ref: https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/Principal-Component-Analysis/Principal-component-analysis-in-R/index.html
# http://strata.uga.edu/8370/lecturenotes/principalComponents.html
metabo[ , c(3:ncol(metabo))] <- apply(metabo[ , c(3:ncol(metabo))], 2, as.numeric)

metabo.pca = prcomp(x = metabo[ , c(3:ncol(metabo))], center = TRUE, scale = TRUE)

# How much variance each principal component explains in the data
var_explained <- data.frame("Var" = 100*metabo.pca$sdev^2/sum(metabo.pca$sdev^2),
                            "PCA" = paste0("PCA", "_", c(1:31))) 
# Scree plot
(scree <- fviz_screeplot(metabo.pca, addlabels = TRUE,
                         barfill = "gray", barcolor = "white",
                         ylim = c(0, 30)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "cm")))

# PCA scatter plot
metabo$group <- factor(metabo$group, levels = c("healthy", "CKD", "PD", "HD", "KTx"))
(PCA_plot <- fviz_pca_ind(metabo.pca,
                          pointsize = 3,
                          label = "none", # hide individual labels
                          habillage = metabo$group, # color by groups
                          palette = "Set2",
                          ellipse.type = "convex",
                          addEllipses = TRUE) +
    scale_shape_manual(values=c(20, 20, 20, 20, 20)) +
    scale_size(1) + 
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          axis.text = element_text(size = 12),
          plot.margin = unit(c(1,1,1,1), "cm")))

# Biplot
(biplot <- fviz_pca_var(metabo.pca, col.var="contrib",
                        labelsize = 2,
                        gradient.cols = c("#00AFBB", "#E7B800", "red"),
                        repel = TRUE, max.overlaps = 20) + 
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          axis.text = element_text(size = 12),
          plot.margin = unit(c(1,1,1,1), "cm")))

# Loading plot
# Contributions of variables to PC1
contr1 <- fviz_contrib(metabo.pca, choice = "var", 
                       axes = 1, top = 10, fill = "gray", color = "white") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        title = element_text(size = 16),
        axis.text = element_text(size = 11),
        plot.margin = unit(c(1,1,1,1), "cm"))
# Contributions of variables to PC2
contr2 <- fviz_contrib(metabo.pca, choice = "var", 
                       axes = 2, top = 10, fill = "gray", color = "white") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        title = element_text(size = 16),
        axis.text = element_text(size = 11),
        plot.margin = unit(c(1,1,1,1), "cm"))

pdf("PCA_analysis.pdf", width = 8, height = 8)
PCA_plot
scree
biplot
contr1
contr2
dev.off()
