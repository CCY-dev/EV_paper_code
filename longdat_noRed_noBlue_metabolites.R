rm(list = ls())
setwd("EV_CKD/20220215_version/Analysis/Longitudinal/Metabolite_as_feature_EV_as_metavariable")
library(metadeconfoundR)
library(readxl)
library(tidyverse)
library(ggrepel)

# Read in the table
data <- read_xlsx("~/Documents/Lab/EV_CKD/20220215_version/Data/Database_CKD-EVs_polished_noRed_noBlue.xlsx",
                  col_names = T)


# Extract only the paired ID data for longitudinal analysis
data_KTx <- data %>% filter(group == "KTx")
data_nonKTx <- data %>% filter(group != "KTx")
paired <- intersect(data_KTx$ID, data_nonKTx$ID)

data <- data %>%
  filter(ID %in% paired)

# Replace "NA"s with NAs
data[data== "NA"] <- NA

# Add new column for test_var
# not KTx = V1 (baseline)
# KTx = V2
data <- data %>%
  mutate(Status = ifelse(data$group == "KTx",
                         yes = "V2", no = "V1"), .before = 2)

colnames(data)[1] <- "Individual"

# Run longdat
library(LongDat)
a = longdat_disc(input = data,
                 test_var = "Status",
                 data_type = "measurement",
                 variable_col = 62,
                 fac_var <- c(1:4, 6:25, 27:37, 48),
                 not_used <- c(3))

result_table = a[[1]]
confound_table = a[[2]]

#write.table(result_table, "20220207_longdatdisc_noRed_noBlue_metabolite_as_feature_EV_as_metavariable_result_table.txt", quote = F, sep = "\t", col.names = T, row.names = F)
#write.table(confound_table, "20220207_longdatdisc_noRed_noBlue_metabolite_as_feature_EV_as_metavariable_confound_table.txt", quote = F, sep = "\t", col.names = T, row.names = F)

result_table <- read.table("20220207_longdatdisc_noRed_noBlue_metabolite_as_feature_EV_as_metavariable_result_table.txt", sep = "\t", header = T)
confound_table <- read.table("20220207_longdatdisc_noRed_noBlue_metabolite_as_feature_EV_as_metavariable_confound_table.txt", sep = "\t", header = T)

result_table$Feature <- c("2-aminophenol", "3-hydroxyanthranilate",
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


# Plotting cuneiform plot
pdf("20220215_longdatdisc_noRed_noBlue_metabolite_as_feature_cuneiform.pdf", width = 8, height = 6)
x_axis_order = NULL
confound_panel = TRUE
pos_color = "gold"
neg_color = "blue3"
panel_width = 4
title = "LongDat result cuneiform plot"
title_size = 16
confound_text_size = 4
x_label_size = 10
y_label_size = 10
legend_title_size = 12
legend_text_size = 10

sig_result <- result_table %>%
  dplyr::filter(Signal != "NS")

if (nrow(sig_result) == 0) {
  stop('All results are insignificant! There is nothing to plot.')}

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

# Plotting
g1 <- ggplot2::ggplot(All_long, aes(x = Effect_name, y = reorder(Feature, EffectSize))) +
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
  labs(fill = "Effect size", alpha = "Significance",
       caption = "V1 = before KTx, V2 = after KTx") +
  theme_light() +
  theme(title = element_text(size = title_size),
        axis.text.y = element_text(size = y_label_size),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size = x_label_size),
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
          axis.text.x=element_text(size = x_label_size),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.x = element_blank())

  final_plot <- (g1|g2) + patchwork::plot_layout(guides = "collect",
                                                 widths = c(panel_width, 1))
  final_plot
dev.off()

# Plotting volcano plot
pdf("20220215_longdatdisc_noRed_noBlue_metabolite_as_feature_volcano.pdf", width = 8, height = 6)
result_extract <- result_table[ , c(1, 4, 6, 7)]
result_extract$logq <- -(log(result_extract$Null_time_model_q, 10))
result_extract$Color <- ifelse(result_extract$Signal == "NS", "Insignificant", "Significant")

ggplot(result_extract, aes(EffectSize_V1_V2, logq)) +
  geom_point(aes(color = Color), alpha = 0.8, size = 2) +
  geom_text_repel(data = subset(result_extract, Signal != "NS"), aes(label = Feature), size = 3, min.segment.length = 0.1) +
  scale_color_manual(values = c("Insignificant" = "darkgrey", "Significant" = "red")) +
  theme_light() +
  labs(x = "Effect size",
       y = "-log(q)") +
  geom_hline(yintercept = 2, linetype="dashed", size = 0.5) +
  geom_vline(xintercept = 0, linetype="dashed", size = 0.5)
dev.off()

