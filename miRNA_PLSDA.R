rm(list = ls())
setwd("EV_CKD/20220406_version/Analysis/All_noRed_noBlue/miRNA_plsda/")
library(tidyverse)
library(caret)
library(readxl)
library(klaR)
library(ggConvexHull)
library(paletteer)

raw <- read.table("Lab/EV_CKD/20220406_version/data/Database_CKD-EVs_polished_noRed_noBlue_miRNA_total_raw.txt",
                  header = T, sep = "\t")
mirna <- raw[ , c(1, 2, 92:ncol(raw))]

# Replace "NA"s with NAs and remove the row with NAs
mirna[mirna == "NA"] <- NA
sum(is.na(mirna))

sort(rowSums(mirna[3:ncol(mirna)]))

# Remove the samples with 0 depth
rownames(mirna) <- paste0(mirna$group, "_", mirna$ID)
values <- mirna[ , -c(1:2)]
values_nonzero <- values[rowSums(values) > 0, ]
  
values_nonzero <- values_nonzero %>% 
  rownames_to_column("IDgroup")
meta <- mirna %>% 
  rownames_to_column("IDgroup") %>% 
  dplyr::select(1:3)

new <- inner_join(meta, values_nonzero, by = "IDgroup")

# Log transform the values
values_new <- new[ , c(4:ncol(new))]
values_new <- apply(values_new, 2, as.numeric)
values_log <- log(values_new + 1)

### Do PLSDA
response = new$group
pls_obj <- plsda(x = values_log, y = as.factor(response), ncomp = 70, probMethod = "Bayes")

# make predictions on the same data
#predictions <- predict(pls_obj, values)
# summarize accuracy
#View(table(predictions, response))

dfscores <- data.frame(Comp1 = pls_obj$scores[ , 1],
                       Comp2 = pls_obj$scores[ , 2],
                       Group = new$group, 
                       ID = new$ID)

dfscores$Group <- factor(dfscores$Group, levels = c("healthy", "CKD", "PD", "HD", "KTx"))

# Plotting
p <- ggplot(dfscores, aes(x=Comp1, y=Comp2)) +
  geom_point(alpha = 0.7, aes(color = Group)) +
  geom_convexhull(aes(fill = Group, color = Group),
                  alpha = 0.2) +
  scale_colour_paletteer_d("RColorBrewer::Set2") + 
  scale_fill_paletteer_d("RColorBrewer::Set2") +
  theme_light()

p

pdf("miRNA_plsda_LogTransformedRawCount_analysis.pdf", width = 8, height = 6)
p
dev.off()





  
