rm(list = ls())
setwd("EV_CKD/20220414_version/analysis/Cross-sectional/")
library(LongDat)
library(readxl)
library(tidyverse)

# Read in the table
data <- read.table("EV_CKD/20220401_version/data/Database_CKD-EVs_polished_noRed_noBlue_miRNA_total_raw.txt",
                   header = T, sep = "\t")
depth_raw <- read.table("EV_CKD/20220302_version/miRNA_depth/Raw_mirna_depth.txt",
                        header = T, sep = "\t")

# Exclude EVs
data <- data[ , -c(48:60)]

# Replace "NA"s with NAs
data[data== "NA"] <- NA

colnames(data)[1] <- "Individual"

value <- data[ , 79:ncol(data)]
value <- apply(value, 2, as.numeric)
val_nonzero <- value[ , colSums(value) > 0]

data_nonzero <- cbind(data[ , 1:78], val_nonzero)

# Add depth info
colnames(depth_raw)[2] <- "Individual"
data_nonzero_final <- inner_join(depth_raw, data_nonzero, by = c("Individual", "group"))
data_nonzero_final$Sample <- NULL

######### Use LongDat scripts to run cross-sectional analysis
library(LongDat)
library(lme4)
library(tidyverse)
library(reshape2)
library(glmmTMB)
library(emmeans)
library(bestNormalize)
library(MASS)

set.seed(100)
input = data_nonzero_final
data_type = "count"
test_var = "group"
variable_col = 80
fac_var = c(1, 2, 4, 6:37, 48)
not_used =  NULL
adjustMethod = "fdr"
model_q = 0.1
posthoc_q = 0.05
theta_cutoff = 2^20
nonzero_count_cutoff1 = 9
nonzero_count_cutoff2 = 5
verbose = T

source("Lab/R_project/longdat/R/data_preprocess.R")
source("Lab/R_project/longdat/R/fix_name_fun.R")
source("Lab/R_project/longdat/R/factor_p_cal.R")
source("Lab/R_project/longdat/R/NuModelTest_disc.R")
source("Lab/R_project/longdat/R/ConModelTest_disc.R")
source("Lab/R_project/longdat/R/unlist_table.R")
source("Lab/R_project/longdat/R/cliff_cal.R")
source("Lab/R_project/longdat/R/wilcox_posthoc.R")
source("Lab/R_project/longdat/R/rm_sparse_disc.R")
source("Lab/R_project/longdat/R/final_result_summarize_disc.R")
source("Lab/R_project/longdat/R/random_neg_ctrl_disc.R")

############## Data preprocessing #################
if (verbose == TRUE) {print("Start data preprocessing.")}
preprocess_lists <- data_preprocess(input, test_var, variable_col,
                                    fac_var, not_used)
mean_abundance <- preprocess_lists[[1]]
prevalence <- preprocess_lists[[2]]
variables_original <- preprocess_lists[[3]]
melt_data <- as.data.frame((preprocess_lists[[4]]))
variables <- preprocess_lists[[5]]
factor_columns <- preprocess_lists[[6]]
factors <- preprocess_lists[[7]]
data <- preprocess_lists[[8]]
values <- preprocess_lists[[9]]
N <- length (variables)
if (verbose == TRUE) {print("Finish data preprocessing.")}

########## Calculate the p values for every factor
#          (used for selecting factors later)
if (variable_col-1-2-length(not_used) > 0) {
  if (verbose == TRUE) {print("Start selecting potential confounders.")}
  factor_p_lists <- suppressWarnings(factor_p_cal(melt_data, variables,
                                                  factor_columns, factors,
                                                  data, N, verbose))
  Ps <- as.data.frame(factor_p_lists[[1]])
  Ps_effectsize <- as.data.frame(factor_p_lists[[2]])
  sel_fac <- factor_p_lists[[3]]
  if (verbose == TRUE) {print("Finished selecting potential confounders.")}
}

################## Null Model Test and Post-Hoc Test #################
if (verbose == TRUE) {print("Start null model test and post-hoc test.")}
NuModel_lists <- NuModelTest_disc(N, data_type, test_var, melt_data,
                                  variables, verbose)
Ps_null_model <- as.data.frame(NuModel_lists[[1]])
Ps_poho_fdr <- as.data.frame(NuModel_lists[[2]])
case_pairs_name <- NuModel_lists[[3]]
if (verbose == TRUE) {print("Finish null model test and post-hoc test.")}

################## Confounding model test ###############
if (variable_col-1-2-length(not_used) > 0) {
  if (verbose == TRUE) {print("Start confounding model test.")}
  ConModel_lists <- ConModelTest_disc(N, variables, melt_data, sel_fac,
                                      data_type, test_var, verbose)
  Ps_conf_model <- ConModel_lists[[1]]
  Ps_inv_conf_model <-ConModel_lists[[2]]
  if (verbose == TRUE) {print("Finish confounding model test.")}
}

####### Unlist the Ps_conf_model and Ps_inv_conf_model ########
suppressWarnings(
  if (variable_col-1-2-length(not_used) > 0) {
    if (verbose == TRUE) {print(
      "Start unlisting tables from confounding model result.")}
    Ps_conf_model_unlist <- unlist_table(Ps_conf_model, N, variables)
    Ps_conf_inv_model_unlist <- unlist_table(Ps_inv_conf_model, N, variables)
    if (verbose == TRUE) {print(
      "Finish unlisting tables from confounding model result.")}
  })

############## Calculate cliff's delta for all variables #################
# Here the code is changed to fit cross-sectional analysis
library(orddom)
case_pairs <- combn(sort(unique(melt_data[ , test_var])), m = 2)
delta <- data.frame(matrix(nrow = length(row.names(Ps_poho_fdr)),
                           ncol = ncol(case_pairs)))
case_pairs_name <- c()
for (i in seq_len(length(row.names(Ps_poho_fdr)))) {
  # loop through all variables
  if (verbose == TRUE) {print(i)}
  bVariable <- variables[i]
  subdata_pre <- subset(melt_data, variable == bVariable)
  counts <- subdata_pre %>% dplyr::count(.data$Individual)
  
  subdata2 <- subdata_pre
  
  for (k in seq_len(ncol(case_pairs))) { # loop through each case pair
    sub3 <- subdata2[subdata2[ , test_var] == case_pairs[1,k], ] %>%
      dplyr::arrange(Individual)
    sub4 <- subdata2[subdata2[ , test_var] == case_pairs[2,k], ] %>%
      dplyr::arrange(Individual)
    
    # Here it's unpaired cliff's delta!
    d <- as.numeric(dmes(x = sub3$value, y = sub4$value)$dc)
    delta[i, k] <- d
    
    name <- paste(case_pairs[1,k], sep = "_", case_pairs[2,k])
    case_pairs_name <- c(case_pairs_name, name)
  }
}
row.names(delta) <- row.names(Ps_poho_fdr)
colnames(delta) <- paste("effect_size", sep = "_", unique(case_pairs_name))

##### Reset the variable names and confounder names to oringinal ####
rownames(Ps_null_model) <- variables_original
rownames(Ps_poho_fdr) <- variables_original
rownames(delta) <- variables_original
variables <- variables_original

if (data_type == "count" & variable_col-1-2-length(not_used) > 0) {
  names(sel_fac) <- variables_original
  rownames(Ps_conf_model_unlist) <- variables_original
  rownames(Ps_conf_inv_model_unlist) <- variables_original
}


######################### Remove the excluded ones #########################
if (data_type == "count") {
  if (verbose == TRUE) {print(
    "Start removing the dependent variables to be exlcuded.")}
  rm_sparse_lists <- rm_sparse_disc(values, data, nonzero_count_cutoff1,
                                    nonzero_count_cutoff2, theta_cutoff,
                                    Ps_null_model, prevalence,
                                    absolute_sparsity, mean_abundance,
                                    Ps_poho_fdr, delta)
  prevalence <- rm_sparse_lists[[1]]
  absolute_sparsity <- rm_sparse_lists[[2]]
  mean_abundance <- rm_sparse_lists[[3]]
  Ps_null_model <-  rm_sparse_lists[[4]]
  Ps_poho_fdr <- rm_sparse_lists[[5]]
  delta <- rm_sparse_lists[[6]]
  variables <- rm_sparse_lists[[7]]
  bac_include <-  rm_sparse_lists[[8]]
  N <- length(variables)
  if (variable_col-1-2-length(not_used) > 0) {
    sel_fac <-  sel_fac[match(bac_include, table = names(sel_fac))]
    Ps_conf_model_unlist <- Ps_conf_model_unlist %>%
      tibble::rownames_to_column() %>%
      dplyr::filter(.data$rowname %in% bac_include) %>%
      tibble::column_to_rownames()
    Ps_conf_inv_model_unlist <- Ps_conf_inv_model_unlist %>%
      tibble::rownames_to_column() %>%
      dplyr::filter(.data$rowname %in% bac_include) %>%
      tibble::column_to_rownames()
  }
  if (verbose == TRUE) {print(
    "Finish removing the dependent variables to be exlcuded.")}
}

################## FDR correction on Ps_null_model (model test) ############
# Do FDR correction (correct for the number of features)
if (verbose == TRUE) {print(
  "Start multiple test correction on null model test p values.")}
adjust_fun <- function(x) p.adjust(p = x, method = adjustMethod)
suppressWarnings(
  Ps_null_model_fdr <- apply(X = Ps_null_model, MARGIN = 2,
                             FUN = adjust_fun))
if (class(Ps_null_model_fdr)[1] != "numeric") {
  Ps_null_model_fdr <- as.data.frame(Ps_null_model_fdr[ , 1])
} else {
  Ps_null_model_fdr <- as.data.frame(Ps_null_model_fdr[1])
  rownames(Ps_null_model_fdr) <- rownames(Ps_null_model)
}
#Reorder the rows of p_lm_fdr to make it follow the order of variables
Ps_null_model_fdr <-
  as.data.frame(Ps_null_model_fdr[match(variables,
                                        rownames(Ps_null_model_fdr)), ])
rownames(Ps_null_model_fdr) <- variables
colnames(Ps_null_model_fdr) <- c("Ps_null_model_fdr")
if (verbose == TRUE) {print(
  "Finish multiple test correction on null model test p values.")}


############## Generate result table as output #################
if (verbose == T) {print("Start generating result tables.")}
results <- final_result_summarize_disc(variable_col, N, Ps_conf_inv_model_unlist, variables, sel_fac, Ps_conf_model_unlist,
                                       model_q, posthoc_q, Ps_null_model_fdr, Ps_null_model, delta, case_pairs, prevalence,
                                       mean_abundance, Ps_poho_fdr, not_used, Ps_effectsize, case_pairs_name, data_type,
                                       false_pos_count = 0, p_wilcox_final)

result_table = results[[2]]
confound_table = results[[1]]

write.table(result_table, "20220420_longdat_code_cross_sectional_raw_miRNA_as_feature_depth_as_confounder_result_table.txt", quote = F, sep = "\t", col.names = T, row.names = F)
write.table(confound_table, "20220420_longdat_code_cross_sectional_raw_miRNA_as_feature_depth_as_confounder_confound_table.txt", quote = F, sep = "\t", col.names = T, row.names = F)

#save(results, file = "20220420_longdat_code_cross_sectional_raw_miRNA_as_feature_depth_as_confounder.RData")
#load("20220420_longdat_code_cross_sectional_raw_miRNA_as_feature_depth_as_confounder.RData")

result_table <- read.table("20220420_longdat_code_cross_sectional_raw_miRNA_as_feature_depth_as_confounder_result_table.txt", sep = "\t", header = T)
confound_table <- read.table("20220420_longdat_code_cross_sectional_raw_miRNA_as_feature_depth_as_confounder_confound_table.txt", sep = "\t", header = T)


### Prepare for plotting
# Remove those group pairs without significance
result_table_new <- result_table[ , c(1, 4:24)]
nonsig_ones <- which(colnames(result_table_new) %in% c("EffectSize_CKD_KTx",
                                                       "EffectSize_CKD_PD", "EffectSize_HD_PD",
                                                       "EffectSize_healthy_KTx",
                                                       "Effect_CKD_KTx",
                                                       "Effect_CKD_PD", "Effect_HD_PD",
                                                       "Effect_healthy_KTx"))
result_table_new <- result_table_new[ , -nonsig_ones]

# Reverse the effect sizes to fit the group order of previous plots
# Here the Effect columns' direction may not be correct (cuz some effect sizes are reversed)
# But Effect columns are just for judging sig/non-sig
final_result <- data.frame(Feature = result_table_new$Feature,
                           Signal = result_table_new$Signal,
                           
                           EffectSize_HD_CKD = result_table_new$EffectSize_CKD_HD,
                           EffectSize_CKD_Healthy = result_table_new$EffectSize_CKD_healthy * (-1),
                           EffectSize_HD_Healthy = result_table_new$EffectSize_HD_healthy * (-1),
                           EffectSize_KTx_HD = result_table_new$EffectSize_HD_KTx,
                           EffectSize_PD_Healthy = result_table_new$EffectSize_healthy_PD,
                           EffectSize_KTx_PD = result_table_new$EffectSize_KTx_PD * (-1),
                           
                           Effect_HD_CKD = result_table_new$Effect_CKD_HD,
                           Effect_CKD_Healthy = result_table_new$Effect_CKD_healthy,
                           Effect_HD_Healthy = result_table_new$Effect_HD_healthy,
                           Effect_KTx_HD = result_table_new$Effect_HD_KTx,
                           Effect_PD_Healthy = result_table_new$Effect_healthy_PD,
                           Effect_KTx_PD  = result_table_new$Effect_KTx_PD
                           )




## Cuneiform plot
x_axis_order = c("Effect_CKD_Healthy", "Effect_PD_Healthy", "Effect_HD_Healthy",
                 "Effect_HD_CKD", "Effect_KTx_PD", "Effect_KTx_HD")
confound_panel = TRUE
pos_color = "gold"
neg_color = "blue3"
panel_width = 4
title = "miRNA result cuneiform plot"
title_size = 16
confound_text_size = 4
x_label_size = 10
y_label_size = 10
legend_title_size = 12
legend_text_size = 10

sig_result <- final_result %>%
  dplyr::filter(Signal != "NS")


# Select the required columns
sig_wide <- sig_result %>%
  dplyr::select(c(Feature, Signal,
                  stringr::str_which(string = colnames(sig_result),
                                     pattern = "Effect")))

# Divide them into 2 tables: Effect and EffectSize
Effect_wide <- sig_wide %>%
  dplyr::select(stringr::str_which(string = colnames(sig_wide),
                                   pattern = "EffectSize",
                                   negate = TRUE))
EffectSize_wide <- sig_wide %>%
  dplyr::select(c(Feature, Signal,
                  stringr::str_which(string = colnames(sig_wide),
                                     pattern = "EffectSize")))

# Transform into long format
Effect_long <- Effect_wide %>%
  tidyr::pivot_longer(stringr::str_which(string = colnames(Effect_wide),
                                         pattern = "Effect"),
                      names_to = "Effect_name", values_to = "Effect")
EffectSize_long <- EffectSize_wide %>%
  tidyr::pivot_longer(stringr::str_which(string = colnames(EffectSize_wide),
                                         pattern = "EffectSize"),
                      names_to = "EffectSize_name", values_to = "EffectSize")

# Merge the two long tables
All_long <- cbind(Effect_long, EffectSize_long[ , c(3:4)])

# Define plotting parameters
All_long$Alpha <- ifelse(All_long$Effect != "NS",
                         yes = "Significant", no = "Insignificant")
All_long$Shape <- ifelse(All_long$EffectSize > 0,
                         yes = "24",
                         no = ifelse(All_long$EffectSize < 0,
                                     yes = "25", no = "1"))

# Reorder plotting sequence
if (is.null(x_axis_order)) {
  All_long$Effect_name <-  All_long$Effect_name
} else {
  All_long$Effect_name <- factor(All_long$Effect_name,
                                 levels = x_axis_order)
}

# Specify the levels of alpha and shape
All_long$Alpha <- factor(All_long$Alpha, levels = c("Significant",
                                                    "Insignificant"))
All_long$Shape <- factor(All_long$Shape, levels = c("1", "24", "25"))

write.table(All_long, "20220420_LongDat_code_cross_sectional_raw_minRNA_as_feature_depth_as_confounder_for_cuneiform_plotting.txt", 
            sep = "\t", quote = F, row.names = F)

# Plotting
#pdf("20220420_longdat_code_cross_sectional_raw_miRNA_as_feature_depth_as_confounder_cuneiform.pdf", 
#    width = 8, height = 10)
g1 <- ggplot2::ggplot(All_long, aes(x = Effect_name, y = Feature)) +
  geom_point(aes(shape = Shape, fill = EffectSize,
                 alpha = Alpha), size = 3.5) +
  scale_shape_manual(values = c(1, 24, 25),
                     breaks = c("1", "24", "25"),
                     labels = c("No change", "Enriched", "Decreased"),
                     name = "Effect", drop = FALSE) +
  scale_fill_gradient2(midpoint = 0, low = neg_color, mid = "white",
                       high = pos_color, n.breaks = 8,
                       limits = c(-1, 1) * max(abs(All_long$EffectSize))) +
  scale_alpha_manual(breaks = c("Insignificant", "Significant"),
                     values=c(0.4, 1), drop = FALSE) +
  ggtitle(title) +
  labs(fill = "Effect size", alpha = "Significance") +
  scale_x_discrete(labels=c("Effect_CKD_Healthy" = "CKD_Healthy", "Effect_PD_Healthy" = "PD_Healthy",
                            "Effect_HD_Healthy" = "HD_Healthy", "Effect_HD_CKD" = "HD_CKD",
                            "Effect_KTx_PD" = "KTx_PD", "Effect_KTx_HD" = "KTx_HD")) +
  theme_light() +
  theme(title = element_text(size = title_size),
        axis.text.y = element_text(size = y_label_size),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size = x_label_size, angle = 90),
        axis.title.y=element_blank(),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size)) +
  guides(fill = guide_colorbar(raster = FALSE, nbin = 30))

g2 <- ggplot2::ggplot(Effect_wide, aes(x = "Confounding status",
                                       y = Feature)) +
  geom_text(aes(label = Signal),
            size = confound_text_size, color = "gray30") +
  theme_light() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size = x_label_size, angle = 90),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.x = element_blank())

final_plot <- (g1|g2) + patchwork::plot_layout(guides = "collect",
                                               widths = c(panel_width, 1))
final_plot
dev.off()

#### Heatmap
All_long$mark <- ifelse(All_long$Effect == "NS", 
                        yes = " ", no = "*")

# Order by name but not clustered automatically
All_long$Effect_name <- factor(All_long$Effect_name, levels = x_axis_order)

# Make marks for heatmap
All_long <- All_long %>% 
  mutate(Mark_new = case_when(
    Effect == "NS" ~ " ",
    Effect != "NS" & Signal == "OK_nc" ~ "o",
    Effect != "NS" & Signal %in% c("C", "AD") ~ "x"
  ))

#pdf("20220420_longdat_code_cross_sectional_raw_miRNA_as_feature_depth_as_confounder_unclustered_heatmap.pdf", width = 6, height = 8)
ggplot(All_long, aes(Effect_name, Feature, fill = EffectSize)) +
  geom_tile(color = "gray88") +
  geom_text(size = 3, aes(label = Mark_new), face = "bold", color = "white") +
  labs(title = " MiRNA result heatmap", 
       subtitle = "o: q < 0.1 and no confounder 
x: q < 0.1 but confounded or ambiguously deconfounded",
       x = "") +
  theme_light() +
  scale_x_discrete(labels=c("Effect_CKD_Healthy" = "CKD_Healthy", "Effect_PD_Healthy" = "PD_Healthy",
                            "Effect_HD_Healthy" = "HD_Healthy", "Effect_HD_CKD" = "HD_CKD",
                            "Effect_KTx_PD" = "KTx_PD", "Effect_KTx_HD" = "KTx_HD")) +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.y = element_text(size = 8)) +                 
  guides(fill = guide_colourbar(title = "Cliff's delta",
                                raster = FALSE, nbin = 25,
                                title.position = "top",  
                                barwidth = 0.5,                               
                                title.hjust = 0.5)) +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "gold") + 
  ggpubr::rotate_x_text(angle = 90)

dev.off()

# Clustered automatically by hclust
# Ref: https://stackoverflow.com/questions/25528059/cluster-data-in-heat-map-in-r-ggplot
final_result2 <- final_result[ , c(1, 3:8)] %>% 
  column_to_rownames("Feature")

pd <- as.data.frame(scale((final_result2)))
ord <- hclust(dist(pd, method = "euclidean"), method = "ward.D2")$order
ord

pd2 <- as.data.frame(scale(t(final_result2)))
ord2 <- hclust(dist(pd2, method = "euclidean"), method = "ward.D")$order
ord2

All_long2 <- All_long
All_long2$Feature <- factor(All_long2$Feature , levels = rownames(final_result2)[ord])
#all_long2$Group <- factor(all_long2$Group, levels = colnames(all_logfold2)[ord2])

# Specify group order manually
#all_long2$Group <- factor(all_long2$Group, levels = c("HD_Healthy", "PD_Healthy", "CKD_Healthy", 
#                                                      "KTx_PD", "KTx_HD"))

# Make marks for heatmap
All_long2 <- All_long2 %>% 
  mutate(Mark_new = case_when(
    Effect == "NS" ~ " ",
    Effect != "NS" & Signal == "OK_nc" ~ "o",
    Effect != "NS" & Signal %in% c("C", "AD") ~ "x"
  ))

#pdf("20220420_longdat_code_cross_sectional_raw_miRNA_as_feature_depth_as_confounder_clustered_heatmap.pdf", width = 6, height = 8)
ggplot(All_long2, aes(Effect_name, Feature, fill = EffectSize)) +
  geom_tile(color = "gray88") +
  geom_text(size = 3, aes(label = Mark_new), color = "white", face = "bold") +
  labs(title = " MiRNA result heatmap", 
       subtitle = "o: q < 0.1 and no confounder 
x: q < 0.1 but confounded or ambiguously deconfounded",
       x = "") +
  theme_light() +
  scale_x_discrete(labels=c("Effect_CKD_Healthy" = "CKD_Healthy", "Effect_PD_Healthy" = "PD_Healthy",
                            "Effect_HD_Healthy" = "HD_Healthy", "Effect_HD_CKD" = "HD_CKD",
                            "Effect_KTx_PD" = "KTx_PD", "Effect_KTx_HD" = "KTx_HD")) +
  theme(plot.title = element_text(size = 18), 
        plot.subtitle = element_text(size = 10), 
        legend.position="right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        axis.text.y = element_text(size = 8)) +                 
  guides(fill = guide_colourbar(title = "Cliff's delta",
                                raster = FALSE, nbin = 25,
                                title.position = "top",  
                                barwidth = 0.5,                               
                                title.hjust = 0.5)) +
  scale_fill_gradient2(low = "blue3", mid = "white", high = "gold") + 
  ggpubr::rotate_x_text(angle = 90)


dev.off()

## Compare the result of longdat and deseq2
library(ggvenn)
deseq_re <- read.table("Lab/EV_CKD/20220406_version/analysis/All_noRed_noBlue/Total_raw_miRNA_GroupToGroup_comparison/20220406_Group_to_group_miRNA_total_raw_Deseq2.txt",
                       sep = "\t", header = T)
longdat_re <- All_long

deseq_sig <- unique(deseq_re$miRNA[deseq_re$q < 0.1])
longdat_sig <- unique(longdat_re$Feature[longdat_re$Effect != "NS"])

for_venn <- list(LongDat = longdat_sig, 
                 Deseq2 = deseq_sig)
# Venn diagram
# Only count the miRNAs, not considering the groups
g3 <- ggvenn(for_venn, stroke_size = 0.2, set_name_size = 4,
             text_size = 4, 
             fill_color = c("#FAC866", "#99C1FA"), auto_scale = T
) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title = element_text(size = 6)) +
  labs(title = "LongDat and Deseq2 comparison: significant miRNA")
g3
ggsave(filename = "LongDat_Deseq2_Venn_miRNA.pdf", device = "pdf",
       width = 7, height = 7)

# Venn diagram
# Count both the miRNAs and the groups
deseq_re <- deseq_re %>% 
  mutate(ID = paste0(miRNA, ",", Group))
longdat_re <- longdat_re %>% 
  mutate(ID = paste0(Feature, ",", 
                     gsub(x = Effect_name, pattern = "Effect_", replacement = "")))

deseq_sig_add_group <- unique(deseq_re$ID[deseq_re$q < 0.1])
longdat_sig_add_group  <- unique(longdat_re$ID[longdat_re$Effect != "NS"])

for_venn_add_group <- list(LongDat = longdat_sig_add_group, 
                 Deseq2 = deseq_sig_add_group)

g4 <- ggvenn(for_venn_add_group, stroke_size = 0.2, set_name_size = 4,
             text_size = 4, auto_scale = T,
             fill_color = c("#FAC866", "#99C1FA")
) +
  labs(title = "LongDat and Deseq2 comparison: significant miRNA within groups") +
  theme(plot.margin = unit(c(0,0,0, 0), "cm"), 
        axis.title = element_text(size = 6)) 
 
g4
ggsave(filename = "LongDat_Deseq2_Venn_miRNA_added_group.pdf", device = "pdf",
       width = 7, height = 7)

