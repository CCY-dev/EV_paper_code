library(tidyverse)
library(DESeq2)
library("RColorBrewer")

## packages for geneset enrichment 
library(gprofiler2)
library(org.Hs.eg.db)
library(GO.db)



files <- list.files(path = ".", recursive = T, pattern = "featureCounts.txt")

joined <- map(files, ~{(readr::read_delim(.x, delim = "\t", skip = 1 ) %>% 
                          dplyr::select(Geneid, contains("bam")))}) %>% 
  purrr::reduce(., dplyr::left_join, by = 'Geneid')


met <- files %>% 
  as_tibble() %>% 
  tidyr::separate(value, into = c("round", "name"), sep = "/") %>% 
  mutate(sample = stringr::str_extract(name, "^.*_")) %>% 
  mutate(round = stringr::str_extract(round, "_*.$")) %>%  
  mutate(sample = stringr::str_remove(sample, "_")) %>% dplyr::select(-name) %>% 
  mutate(group = case_when(stringr::str_detect(sample, "KG")  ~ "control", # is healthy 
                           stringr::str_detect(sample, "XE")  ~ "veh_control",
                           stringr::str_detect(sample, "HD")  ~ "HD",
                           #stringr::str_detect(sample, "HC")  ~ "HC",
                           stringr::str_detect(sample, "NTX")  ~ "KTx")) %>% 
  mutate(round = case_when(round == "1" ~ "one", 
                           TRUE ~ "two"))


# Remove "Aligned.sortedByCoord.out.bam" from all column names
names(joined) <- str_remove_all(names(joined), "Aligned.sortedByCoord.out.bam") 

joined


colDat <- met %>% 
  filter(!stringr::str_detect(sample, "XE")) %>% 
  column_to_rownames("sample") %>% 
  mutate_if(is.character, as.factor)
cntDat <- joined %>% dplyr::select(-contains("XE")) %>% column_to_rownames("Geneid")

all(rownames(colDat) == colnames(cntDat))



dds <- DESeqDataSetFromMatrix(countData = cntDat, 
                              colData = colDat, 
                              design = ~ group)


## Combat Seq is better # 
library(sva)
# adjusting for the rounds 
adjusted <- ComBat_seq(as.matrix(cntDat), batch=colDat$round, group=colDat$group)

dds2 <- DESeqDataSetFromMatrix(countData = adjusted, 
                               colData = colDat, 
                               design = ~ group)
# pre-filtering 
#keep only rows that have at least 10 reads total.
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]

vsd2 <- vst(dds2, blind=FALSE)



dds2 <- DESeq(dds2)


res <- results(dds2, contrast=c("group","KTx","HD")) # specify what you want to compare this way 
res
## Estimate shrinkage before counting and summarising P values #
# shrinkage of log2fc allows for stringent gene-filtering (reduce false positives)
resLFC <- lfcShrink(dds2, coef="group_HD_vs_control", type="apeglm")

resLFC

# order p val 
resOrdered <- resLFC[order(resLFC$pvalue),]


# how many adjusted p values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

# diff ways to look at it 

res_HD_vs_KTx <- results(dds2, alpha=0.1, contrast=c("group","HD","KTx"))
res_HD_vs_control <- results(dds2, alpha=0.1, contrast=c("group","HD","control"))


res_HD_vs_KTx <- res_HD_vs_KTx %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  mutate(comp = "HD_vs_KTx")


res_HD_vs_control <- res_HD_vs_control %>%  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  mutate(comp = "HD_vs_control")

# GSEA

go_terms <- c("GO:0036120", "GO:0000082", "GO:0048660", "GO:0045765", "GO:0030336", "GO:0001525")

retrieved2 <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=go_terms, columns=c("ENSEMBL", "SYMBOL"))

retrieved2 %>% group_by(GOALL)

ang_genes <- retrieved2 %>% filter(GOALL =="GO:0001525") %>% pull(ENSEMBL)



res_HD_vs_KTx_v <- res_HD_vs_KTx %>% rename("gene" = "ENSEMBL")
res_HD_vs_control_v <- res_HD_vs_control %>% rename("gene" = "ENSEMBL")

met <- colDat %>% rownames_to_column("sam") %>% mutate(sample_name = paste0("X", sam)) %>% as_tibble()


sel <- res_HD_vs_KTx_v %>% filter(log2FoldChange > 0.1 & pvalue < 0.1) %>% distinct(ENSEMBL) %>% 
  full_join(res_HD_vs_control_v %>% filter(log2FoldChange > 0.1& pvalue < 0.1) %>% distinct(ENSEMBL)) %>% 
  distinct()


ok <- sel %>% filter(ENSEMBL %in%retrieved2$ENSEMBL) %>% 
  inner_join(., retrieved2 %>% dplyr::select(ENSEMBL, GOALL, SYMBOL))



ok2 <- ok %>% 
  mutate(GO_Term = case_when(GOALL == "GO:0001525" ~ "Angiogenesis",
                             TRUE ~ "other")) %>% 
  dplyr::select(ENSEMBL, SYMBOL, GO_Term) %>% 
  distinct() %>% 
  group_by(SYMBOL) %>%
  dplyr::slice(if(any(GO_Term == "Angiogenesis")) which.max(GO_Term == "Angiogenesis") else 1) %>%
  ungroup()
# 
# 
# ok %>% distinct() %>% 
#   write_csv("GO_genes.csv")


### Heatmap ######### 

library(tidyHeatmap)
library(ComplexHeatmap)
fexp <- c("#254093", "#2d4dab", "#355bc3", "#3c69da", "#4477f2",  
          "#4d86ff", "#5594ff", "#5da3ff", "#66b2ff", "#6fc1ff",  
          "#ffffff", "#fffdd9", "#fffbaf", "#fff995", "#fff77c", 
          "#fff562", "#fff349", "#fff130", "#ffef17", "#ffde55")

counts(dds2, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="ENSEMBL") %>%
  as_tibble()%>% 
  inner_join(ok2) %>% 
  tidyr::pivot_longer(-c(ENSEMBL,GO_Term, SYMBOL),  names_to = "sample_name", values_to = "ab") %>% 
  inner_join(met) %>% #View()
  mutate(group = as.factor(group),
         GO_Term = as.factor(GO_Term) ) %>% 
  group_by(GO_Term) %>% 
  tidyHeatmap::heatmap(SYMBOL, sample_name, ab, scale = "row",
                       column_names_gp = gpar(fontsize = 5),
                       row_names_gp = gpar(fontsize = 10),
                       palette_value = circlize::colorRamp2(
                         seq(-4, 4, length.out = 11), 
                         RColorBrewer::brewer.pal(11, "BrBG")
                       )
  ) %>% 
  tidyHeatmap::add_tile(group) %>% 
  tidyHeatmap::add_tile(GO_Term)

##### PCA #########
# library(devtools)
# devtools::install_github("cmartin/ggConvexHull")

pcaData_go <- plotPCA(vsd2[ok2 %>% pull(ENSEMBL),], intgroup=c("group", "round"), returnData=TRUE)
percentVar_go <- round(100 * attr(pcaData_go, "percentVar"))

pcaData_go %>%
  mutate(sample_name = paste0("X", name)) %>% 
  inner_join(met, by = "sample_name") %>% 
  ggplot(., aes(PC1, PC2, group = group.1)) +
  geom_jitter(size=4, aes(colour = group.1)) +
  ggConvexHull::geom_convexhull(alpha = 0.3,aes(fill = group.1)) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  ggsci::scale_color_lancet() + 
  ggsci::scale_fill_lancet() + 
  ggpubr::theme_pubr()
