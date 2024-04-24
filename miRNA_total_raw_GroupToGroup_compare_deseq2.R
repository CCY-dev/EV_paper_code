rm(list = ls())
setwd("EV_CKD/20220406_version/analysis/All_noRed_noBlue/Total_raw_miRNA_GroupToGroup_comparison/")
library(tidyverse)
library(readxl)
library(ggpubr)
library(DESeq2)
library(apeglm)
library(ashr)

raw <- read.table("Lab/EV_CKD/20220406_version/data/Database_CKD-EVs_polished_noRed_noBlue_miRNA_total_raw.txt",
                  header = T, sep = "\t")
mirna <- raw[ , c(1, 2, 92:ncol(raw))]
mirna$group <- gsub(pattern = "healthy", replacement = "Healthy", x = mirna$group)

# Replace "NA"s with NAs and remove the row with NAs
mirna[mirna == "NA"] <- NA
sum(is.na(mirna))

rownames(mirna) <- paste0(mirna$group, "_", mirna$ID)

mirna_t <- t(mirna)
counts <- mirna_t[-c(1, 2), ]
meta <- mirna[ , c(1, 2)]

meta$group <- as.factor(meta$group)
mode(counts) <- "numeric"

# Check if colnames(counts) and rownames(meta) have same order
match(rownames(meta), colnames(counts))
all(rownames(meta) %in% colnames(counts))
all(colnames(counts) %in% rownames(meta))

# Add pseudo 1 to each cell in counts or Deseq2 throws error:
# every gene contains at least one zero, cannot compute log geometric means
# https://www.biostars.org/p/440379/
counts_new <- counts + 1

# Create DESeqDataSet for Deseq2
dds <- DESeqDataSetFromMatrix(countData = counts_new,
                              colData = meta,
                              design = ~ group)

# Perform Deseq2
output <- DESeq(dds)

# Get pairwise comparisons
# https://www.biostars.org/p/325009/
CKD_Healthy_shk <- as.data.frame(lfcShrink(output, contrast=c("group", "CKD", "Healthy"), type = "ashr"))
PD_Healthy_shk <-  as.data.frame(lfcShrink(output, contrast=c("group", "PD", "Healthy"), type = "ashr"))
HD_Healthy_shk <-  as.data.frame(lfcShrink(output, contrast=c("group", "HD", "Healthy"), type = "ashr"))
KTx_Healthy_shk <-  as.data.frame(lfcShrink(output, contrast=c("group", "KTx", "Healthy"), type = "ashr"))
PD_CKD_shk <-  as.data.frame(lfcShrink(output, contrast=c("group", "PD", "CKD"), type = "ashr"))
HD_CKD_shk <-  as.data.frame(lfcShrink(output, contrast=c("group", "HD", "CKD"), type = "ashr"))
KTx_CKD_shk <-  as.data.frame(lfcShrink(output, contrast=c("group", "KTx", "CKD"), type = "ashr"))
HD_PD_shk <-  as.data.frame(lfcShrink(output, contrast=c("group", "HD", "PD"), type = "ashr"))
KTx_PD_shk <-  as.data.frame(lfcShrink(output, contrast=c("group", "KTx", "PD"), type = "ashr"))
KTx_HD_shk <-  as.data.frame(lfcShrink(output, contrast=c("group", "KTx", "HD"), type = "ashr"))

#plotMA(lfcShrink(output, contrast=c("group", "HD", "Healthy"), type = "ashr"), ylim=c(-2,2))

# Merge the results
all_logfold <- data.frame(CKD_Healthy = CKD_Healthy_shk$log2FoldChange,
                          PD_Healthy = PD_Healthy_shk$log2FoldChange,
                          HD_Healthy = HD_Healthy_shk$log2FoldChange,
                          KTx_Healthy = KTx_Healthy_shk$log2FoldChange,
                          PD_CKD = PD_CKD_shk$log2FoldChange,
                          HD_CKD = HD_CKD_shk$log2FoldChange,
                          KTx_CKD = KTx_CKD_shk$log2FoldChange,
                          HD_PD = HD_PD_shk$log2FoldChange,
                          KTx_PD = KTx_PD_shk$log2FoldChange,
                          KTx_HD = KTx_HD_shk$log2FoldChange, 
                          row.names = rownames(CKD_Healthy_shk))

all_p <- data.frame(CKD_Healthy = CKD_Healthy_shk$pvalue,
                          PD_Healthy = PD_Healthy_shk$pvalue,
                          HD_Healthy = HD_Healthy_shk$pvalue,
                          KTx_Healthy = KTx_Healthy_shk$pvalue,
                          PD_CKD = PD_CKD_shk$pvalue,
                          HD_CKD = HD_CKD_shk$pvalue,
                          KTx_CKD = KTx_CKD_shk$pvalue,
                          HD_PD = HD_PD_shk$pvalue,
                          KTx_PD = KTx_PD_shk$pvalue,
                          KTx_HD = KTx_HD_shk$pvalue, 
                          row.names = rownames(CKD_Healthy_shk))

# Convert to long format
all_logfold_long <- all_logfold %>% 
  rownames_to_column("miRNA") %>% 
  pivot_longer(2:11, values_to = "logfold", names_to = "Group")

pvalue_long <- all_p %>% 
  rownames_to_column("miRNA") %>% 
  pivot_longer(2:11, values_to = "p", names_to = "Group")
pvalue_long$q <- p.adjust(pvalue_long$p, method = "fdr") # FDR correction
pvalue_long$mark <- ifelse(pvalue_long$q < 0.001, yes = "***", 
                           no = ifelse(pvalue_long$q < 0.01, yes = "**", 
                                       no = ifelse(pvalue_long$q < 0.1,
                                                   yes = "*", no = " ")))
pvalue_long$mark[is.na(pvalue_long$mark)] <- " "
all_long <- inner_join(pvalue_long, all_logfold_long, by = c("miRNA", "Group"))  

# Extract only the sig ones
sig_miRNA <- unique(all_long$miRNA[all_long$q < 0.1])
all_long <- all_long %>% 
  filter(miRNA %in% sig_miRNA)

# Extract only the sig groups
sig_group <- unique(all_long$Group[all_long$q < 0.1])
all_long <- all_long %>% 
  filter(Group %in% sig_group)
write.table(all_long, "20220406_Group_to_group_miRNA_total_raw_Deseq2.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#### Heatmap
# Order by name but not clustered automatically
all_long$Group <- factor(all_long$Group, levels = unique(all_long$Group))

pdf("20220406_Group_to_group_miRNA_total_raw_heatmap_unclustered_Deseq2.pdf", width = 6, height = 6)
ggplot(all_long, aes(Group, miRNA, fill = logfold)) +
  geom_tile(color = "gray88") +
  geom_text(size = 2, aes(label = mark)) +
  labs(title = "MiRNA Group comparison", 
       subtitle = "FDR-corrected p values: ***q< 0.001, **q< 0.01, *q< 0.1") +
  theme_light() +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.y = element_text(size = 8)) +                 
  guides(fill = guide_colourbar(title = "Log2 Fold Change",
                                raster = FALSE, nbin = 25,
                                title.position = "top",  
                                barwidth = 0.5,                               
                                title.hjust = 0.5)) +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "gold",
                       breaks = seq(-2.5, 2, 0.5)) + 
  ggpubr::rotate_x_text(angle = 90)

dev.off()

# Clustered automatically by hclust
# Ref: https://stackoverflow.com/questions/25528059/cluster-data-in-heat-map-in-r-ggplot
all_logfold2 <- all_logfold
all_logfold2 <- all_logfold2[ , c(1, 2, 3, 9, 10)]

pd <- as.data.frame(scale((all_logfold2)))
ord <- hclust(dist(pd, method = "euclidean"), method = "ward.D2")$order
ord

pd2 <- as.data.frame(scale(t(all_logfold2)))
ord2 <- hclust(dist(pd2, method = "euclidean"), method = "ward.D")$order
ord2

all_long2 <- all_long
all_long2$miRNA <- factor(all_long2$miRNA, levels = rownames(all_logfold2)[ord])
#all_long2$Group <- factor(all_long2$Group, levels = colnames(all_logfold2)[ord2])

# Specify group order manually
all_long2$Group <- factor(all_long2$Group, levels = c("HD_Healthy", "PD_Healthy", "CKD_Healthy", 
                                                      "KTx_PD", "KTx_HD"))


pdf("20220406_Group_to_group_miRNA_total_raw_heatmap_clustered_Deseq2.pdf", width = 6, height = 6)
ggplot(all_long2, aes(Group, miRNA, fill = logfold)) +
  geom_tile(color = "gray88") +
  geom_text(size = 2, aes(label = mark)) +
  labs(title = "MiRNA Group comparison", 
       subtitle = "FDR-corrected p values: ***q< 0.001, **q< 0.01, *q< 0.1") +
  theme_light() +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.y = element_text(size = 8)) +                 
  guides(fill = guide_colourbar(title = "Log2 Fold Change",
                                raster = FALSE, nbin = 25,
                                title.position = "top",  
                                barwidth = 0.5,                               
                                title.hjust = 0.5)) +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "gold",
                       breaks = seq(-2.5, 2, 0.5)) + 
  ggpubr::rotate_x_text(angle = 90)

dev.off()


# Check the peculiar hsa.miR.4485 
library(ggpubr)
sub <- cbind(raw[ , 1:2], raw$hsa.miR.4485)
colnames(sub)[3] <- "hsa.miR.4485"

sub$group <- factor(sub$group, levels = c("healthy", "KTx", "CKD", "PD", "HD"))
my_compare <- list(c("KTx", "HD"), c("healthy", "HD"))

ggplot(sub, aes(group, hsa.miR.4485)) +
  geom_boxplot() +
  theme_light() +
  stat_compare_means(method = "kruskal.test", label.y = 25) +
  stat_compare_means(method = "wilcox.test", comparisons = my_compare) +
  geom_jitter(alpha = 0.6, color = "salmon", width = 0.1, height = 0.1)

ggsave("hsa.miR.4485_groups.pdf", device = "pdf", width = 8, height = 6)  
  
  
  


