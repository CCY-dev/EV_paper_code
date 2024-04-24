rm(list = ls())
setwd("EV_CKD/20220302_version/All_noRed_noBlue/Metabolite_GroupToGroup_comparison/")
library(tidyverse)
library(readxl)
library(orddom)
library(ggpubr)

raw <- read_xlsx("../../Data/Database_CKD-EVs_polished_noRed_noBlue.xlsx", 
                 col_names = T)
metabo <- raw[ , c(1, 2, 61:91)]
metabo$group <- gsub(pattern = "healthy", replacement = "Healthy", x = metabo$group)

# Replace "NA"s with NAs and remove the row with NAs
metabo[metabo == "NA"] <- NA
sum(is.na(metabo))
metabo <- metabo %>% 
  drop_na()
metabo_long <- metabo %>% 
  pivot_longer(3:ncol(metabo), names_to = "Metabolite", values_to = "Value")
metabo_long$Value <- as.numeric(metabo_long$Value)
metabolites <- unique(metabo_long$Metabolite)
N = length(metabolites)

# Effect sizes between groups
group_pairs <- data.frame(V1 = c("CKD", "Healthy"),
                           V2 = c("PD", "Healthy"),
                           V3 = c("HD", "Healthy"),
                           V4 = c("KTx", "Healthy"),
                           V5 = c("PD", "CKD"),
                           V6 = c("HD", "CKD"),
                           V7 = c("KTx", "CKD"),
                           V8 = c("HD", "PD"),
                           V9 = c("KTx", "PD"),
                           V10 = c("KTx", "HD"))
delta <- data.frame(matrix(nrow = N,
                           ncol = ncol(group_pairs)))
pvalue <- data.frame(matrix(nrow = N,
                           ncol = ncol(group_pairs)))
for (i in 1:N) {
  # loop through all variables
  bVariable <- metabolites[i]
  subdata <- subset(metabo_long, Metabolite == bVariable)
  
  for (k in 1:ncol(group_pairs)) { # loop through each group pair
    sub3 <- subdata[subdata$group == group_pairs[1,k], ]
    sub4 <- subdata[subdata$group == group_pairs[2,k], ]
    
    # Effect size
    # Note that y is treated while x is control
    # This is a non-paired test, so extract "dc"
    d <- as.numeric(dmes(y = sub3$Value, x = sub4$Value)$dc)
    delta[i, k] <- d
    
    # P-value
    w <- wilcox.test(sub3$Value, sub4$Value, paired = F)
    pvalue[i, k] <- w$p.value
  }
}
row.names(delta) <- metabolites
colnames(delta) <- paste0(group_pairs[1, ], sep = " v.s. ", group_pairs[2, ])
row.names(pvalue) <- metabolites
colnames(pvalue) <- paste0(group_pairs[1, ], sep = " v.s. ", group_pairs[2, ])

# Change names of variables
row.names(delta) <- c("2-aminophenol", "3-hydroxyanthranilate",
                      "3-hydroxykynurenine", "3-indoleacetamide",
                      "3-indoleacetonitril", "3-methoxyanthranilate",
                      "5-hydroxyindoleacetate", "5-hydroxytryptophan",
                      "5-methoxyindoleacetate", "5-methoxytryptamine",
                      "6-hydroxymelatonin", "acetylserotonin",
                      "anthranilate", "formylkynurenine",
                      "indole", "indole-3-carboxaldehyde",
                      "indole-3-carboxylic acid", "indole-3-ethanol",
                      "indole-3-methyl-acetate", "indole-3-propionic",
                      "indole-acetic acid", "indolelactate",
                      "indoxylsulfate", "kynurenic acid",
                      "kynurenine", "methyltryptamine",
                      "quinolinic acid", "serotonin",
                      "tryptamin", "tryptophan",
                      "xanthurenic acid")
row.names(pvalue) <- c("2-aminophenol", "3-hydroxyanthranilate",
                      "3-hydroxykynurenine", "3-indoleacetamide",
                      "3-indoleacetonitril", "3-methoxyanthranilate",
                      "5-hydroxyindoleacetate", "5-hydroxytryptophan",
                      "5-methoxyindoleacetate", "5-methoxytryptamine",
                      "6-hydroxymelatonin", "acetylserotonin",
                      "anthranilate", "formylkynurenine",
                      "indole", "indole-3-carboxaldehyde",
                      "indole-3-carboxylic acid", "indole-3-ethanol",
                      "indole-3-methyl-acetate", "indole-3-propionic",
                      "indole-acetic acid", "indolelactate",
                      "indoxylsulfate", "kynurenic acid",
                      "kynurenine", "methyltryptamine",
                      "quinolinic acid", "serotonin",
                      "tryptamin", "tryptophan",
                      "xanthurenic acid")


# Convert to long format
delta_long <- delta %>% 
  rownames_to_column("Metabolite") %>% 
  pivot_longer(2:11, values_to = "cliffdelta", names_to = "Group")

pvalue_long <- pvalue %>% 
  rownames_to_column("Metabolite") %>% 
  pivot_longer(2:11, values_to = "P", names_to = "Group")
pvalue_long$Q <- p.adjust(pvalue_long$P, method = "fdr")
pvalue_long$mark <- ifelse(pvalue_long$Q < 0.001, yes = "***", 
                           no = ifelse(pvalue_long$Q < 0.01, yes = "**", 
                                       no = ifelse(pvalue_long$Q < 0.1,
                                                   yes = "*", no = " ")))
pvalue_long$mark[is.na(pvalue_long$mark)] <- " "
all_long <- inner_join(pvalue_long, delta_long, by = c("Metabolite", "Group"))  

#### Only plot the first 4 columns
all_long <- all_long %>% 
  filter(Group %in% c("CKD v.s. Healthy", "PD v.s. Healthy", "HD v.s. Healthy",
                      "KTx v.s. Healthy"))
  
#### Heatmap
# Order by name but not clustered automatically
all_long$Group <- factor(all_long$Group, levels = unique(all_long$Group))

pdf("202202302_Group_to_group_metabolite_heatmap_unclustered.pdf", width = 8, height = 8)
ggplot(all_long, aes(Group, Metabolite, fill = cliffdelta)) +
  geom_tile() +
  geom_text(size = 2, aes(label = mark)) +
  labs(title = "Group comparison effect sizes", 
       subtitle = "FDR-corrected p values: ***q< 0.001, **q< 0.01, *q< 0.1") +
  theme_light() +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        #plot.caption = element_text(hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +                 
  guides(fill = guide_colourbar(title = "Cliff's delta",
                                raster = FALSE, nbin = 25,
                                title.position = "top",                  
                                barwidth = 0.5,                               
                                title.hjust = 0.5)) +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "gold") + 
  ggpubr::rotate_x_text(angle = 90)

dev.off()



# Clustered automatically by hclust (only by effect size, don't change the group)
# Ref: https://stackoverflow.com/questions/25528059/cluster-data-in-heat-map-in-r-ggplot
delta2 <- delta

pd <- as.data.frame(scale((delta2)))
ord <- hclust(dist(pd, method = "euclidean"), method = "ward.D2")$order
ord

#pd2 <- as.data.frame(scale(t(delta2)))
#ord2 <- hclust(dist(pd2, method = "euclidean"), method = "ward.D")$order
#ord2

all_long2 <- all_long
all_long2$Metabolite <- factor(all_long2$Metabolite, levels = rownames(delta2)[ord])
#all_long2$Group <- factor(all_long2$Group, levels = colnames(delta2)[ord2])

pdf("20220302_Group_to_group_metabolite_heatmap_clustered.pdf", width = 8, height = 8)
ggplot(all_long2, aes(Group, Metabolite, fill = cliffdelta)) +
  geom_tile() +
  geom_text(size = 2, aes(label = mark)) +
  labs(title = "Group comparison effect sizes", 
       subtitle = "FDR-corrected p values: ***q< 0.001, **q< 0.01, *q< 0.1") +
  theme_light() +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        #plot.caption = element_text(hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +                 
  guides(fill = guide_colourbar(title = "Cliff's delta",
                                raster = FALSE, nbin = 25,
                                title.position = "top",                  
                                barwidth = 0.5,                               
                                title.hjust = 0.5)) +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "gold") + 
  ggpubr::rotate_x_text(angle = 90)

dev.off()



