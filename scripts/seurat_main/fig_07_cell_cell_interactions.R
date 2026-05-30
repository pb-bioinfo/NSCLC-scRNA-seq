library(scales)  # For DiscretePalette function
library(ggplot2)
library(qs)
library(hrbrthemes)
library(patchwork)
library(cowplot)   # For get_legend() and ggdraw()

library(data.table)
library(reshape2)
library(viridis)

library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(ComplexHeatmap)

library(dplyr)
library(Seurat)
library(presto)
library(enrichR)
library(CellChat)

setwd("~/diskD/lung/paneltotal/fig7/")

color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f68231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

cellchat <- qread("cellchat_total.qs")
cellchat_black <- qread("cellchat_black.qs")
cellchat_white <- qread("cellchat_white.qs")


#####
#total incoming and outgoing interaction strength 
extract_group_strength_total <- function(obj, group1, group2) {
  count_mat <- obj@net$count
  weight_mat <- obj@net$weight
  
  # --- Group 1 metrics (e.g., Myeloid Cell Network Profiles) ---
  # Outgoing/Incoming Strength (using weight matrix)
  g1_out_wt <- rowSums(weight_mat[group1, c(1:41), drop = FALSE])
  g1_in_wt  <- colSums(weight_mat[ c(1:41), group1, drop = FALSE])
  
  # BUBBLE SIZE (Using count_mat for the total interactions)
  g1_total_count <- rowSums(count_mat[group1,  c(1:41), drop = FALSE]) + 
    colSums(count_mat[ c(1:41), group1, drop = FALSE])
  
  df1 <- data.frame(labels = group1, Out_Weight = g1_out_wt, In_Weight = g1_in_wt, Total_Count = g1_total_count)
  
  # --- Group 2 metrics (e.g., T Cell Network Profiles) ---
  # Outgoing/Incoming Strength (using weight matrix)
  g2_out_wt <- rowSums(weight_mat[group2,  c(1:41), drop = FALSE])
  g2_in_wt  <- colSums(weight_mat[ c(1:41), group2, drop = FALSE])
  g2_total_count <- rowSums(count_mat[group2,  c(1:41), drop = FALSE]) + 
    colSums(count_mat[ c(1:41), group2, drop = FALSE])
  df2 <- data.frame(labels = group2, Out_Weight = g2_out_wt, In_Weight = g2_in_wt, Total_Count = g2_total_count)
  
  return(rbind(df1, df2))
}



mye_cells <- levels(cellchat@idents[[1]])[21:25] 
t_cells <- levels(cellchat@idents[[1]])[29:39] 



df_black_strength_mye_T <- extract_group_strength_total(cellchat_black, mye_cells, t_cells)
df_white_strength_mye_T <- extract_group_strength_total(cellchat_white, mye_cells, t_cells)


df_diff_strength_mye_T <- merge(df_black_strength_mye_T, df_white_strength_mye_T, by = "labels", suffixes = c("_Black", "_White"))

df_diff_strength_mye_T <- df_diff_strength_mye_T %>%
  mutate(
    Delta_Out_Weight = Out_Weight_Black - Out_Weight_White,
    Delta_In_Weight  = In_Weight_Black - In_Weight_White,
    Avg_Count        = (Total_Count_Black + Total_Count_White) / 2
  )


df_diff_strength_mye_T <- df_diff_strength_mye_T %>%
  mutate(
    Ratio_Out_Weight = log2(Out_Weight_Black / Out_Weight_White),
    Ratio_In_Weight  = log2(In_Weight_Black / In_Weight_White),
  )


f7a <- ggplot(df_diff_strength_mye_T, aes(x = Ratio_Out_Weight, y = Ratio_In_Weight)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  
  geom_point(aes(size = Avg_Count, color = labels), alpha = 0.8) +
  geom_text_repel(aes(label = labels), size = 5, box.padding = 0.5) +
  theme_bw() +
  labs(
    title = "Ratio Signal Strength: Myeloid & T",
    x = "Log2 FC of Outgoing Strength (Black/White)",
    y = "Log2 FC of Incoming Strength (Black/White)",
    size = "Avg Interaction\nCount",
    color = "Cell Type"
  )  +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 12)
  )



ggsave(filename = "f7a.png", plot =f7a, width =7, height = 6.1, units = "in")


#####


extract_group_strength <- function(obj, group1, group2) {
  count_mat <- obj@net$count
  weight_mat <- obj@net$weight
  
  
  g1_out_wt <- rowSums(weight_mat[group1, group2, drop = FALSE])
  g1_in_wt  <- colSums(weight_mat[group2, group1, drop = FALSE])
  
  # BUBBLE SIZE (Using count_mat for the total interactions)
  g1_total_count <- rowSums(count_mat[group1, group2, drop = FALSE]) + 
    colSums(count_mat[group2, group1, drop = FALSE])
  
  df1 <- data.frame(labels = group1, Out_Weight = g1_out_wt, In_Weight = g1_in_wt, Total_Count = g1_total_count)
  
  g2_out_wt <- rowSums(weight_mat[group2, group1, drop = FALSE])
  g2_in_wt  <- colSums(weight_mat[group1, group2, drop = FALSE])
  g2_total_count <- rowSums(count_mat[group2, group1, drop = FALSE]) + 
    colSums(count_mat[group1, group2, drop = FALSE])
  
  df2 <- data.frame(labels = group2, Out_Weight = g2_out_wt, In_Weight = g2_in_wt, Total_Count = g2_total_count)
  
  return(rbind(df1, df2))
}

b_cells <- levels(cellchat@idents[[1]])[1:4] 
t_cells <- levels(cellchat@idents[[1]])[c(29:35,37:39)] 

# 3. Extract the data
df_black_strength_B_T <- extract_group_strength(cellchat_black, b_cells, t_cells)
df_white_strength_B_T <- extract_group_strength(cellchat_white, b_cells, t_cells)


df_diff_strength_B_T <- merge(df_black_strength_B_T, df_white_strength_B_T, by = "labels", suffixes = c("_Black", "_White"))

df_diff_strength_B_T <- df_diff_strength_B_T %>%
  mutate(
    Delta_Out_Weight = Out_Weight_Black - Out_Weight_White,
    Delta_In_Weight  = In_Weight_Black - In_Weight_White,
    Avg_Count        = (Total_Count_Black + Total_Count_White) / 2
  )


df_diff_strength_B_T <- df_diff_strength_B_T %>%
  mutate(
    Ratio_Out_Weight = log2(Out_Weight_Black / Out_Weight_White),
    Ratio_In_Weight  = log2(In_Weight_Black / In_Weight_White),
  )


f7b <- ggplot(df_diff_strength_B_T, aes(x = Ratio_Out_Weight, y = Ratio_In_Weight)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  
  geom_point(aes(size = Avg_Count, color = labels), alpha = 0.8) +
  geom_text_repel(aes(label = labels), size = 5, box.padding = 0.5) +
  
  theme_bw() +
  labs(
    title = "Ratio Signal Strength: B & T Only",
    x = "Log2 FC of Outgoing Strength (Black/White)",
    y = "Log2 FC of Incoming Strength (Black/White)",
    size = "Avg Interaction\nCount",
    color = "Cell Type"
  )  +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 12)
  )

ggsave(filename = "f7b.png", plot = f7b, width = 6, height = 6, units = "in")


#####


mali_cells <- c(levels(cellchat@idents[[1]])[14])
cd8_cells <- levels(cellchat@idents[[1]])[36:39] 

nk_cells <- levels(cellchat@idents[[1]])[26:28] 
cyto_cells <- c(cd8_cells,nk_cells)


# 3. Extract the data
df_black_strength_mali_cyto <- extract_group_strength(cellchat_black, cyto_cells, mali_cells)
df_white_strength_mali_cyto <- extract_group_strength(cellchat_white, cyto_cells, mali_cells)


df_diff_strength_mali_cyto <- merge(df_black_strength_mali_cyto, df_white_strength_mali_cyto, by = "labels", suffixes = c("_Black", "_White"))

df_diff_strength_mali_cyto <- df_diff_strength_mali_cyto %>%
  mutate(
    Delta_Out_Weight = Out_Weight_Black - Out_Weight_White,
    Delta_In_Weight  = In_Weight_Black - In_Weight_White,
    Avg_Count        = (Total_Count_Black + Total_Count_White) / 2
  )


df_diff_strength_mali_cyto <- df_diff_strength_mali_cyto %>%
  mutate(
    Ratio_Out_Weight = log2(Out_Weight_Black / Out_Weight_White),
    Ratio_In_Weight  = log2(In_Weight_Black / In_Weight_White),
  )


f7c <- ggplot(df_diff_strength_mali_cyto, aes(x = Ratio_Out_Weight, y = Ratio_In_Weight)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  
  geom_point(aes(size = Avg_Count, color = labels), alpha = 0.8) +
  geom_text_repel(aes(label = labels), size = 5, box.padding = 0.5) +
  
  theme_bw() +
  labs(
    title = "Ratio Signal Strength: T/NK vs Malignant Cells Only",
    x = "Log2 FC of Outgoing Strength (Black/White)",
    y = "Log2 FC of Incoming Strength (Black/White)",
    size = "Avg Interaction\nCount",
    color = "Cell Type"
  )  +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 12)
  )



ggsave(filename = "f7c.png", plot = f7c, width = 6, height = 6, units = "in")


#####


pathway_name <- "SPP1"

p_black <- netVisual_heatmap(
  cellchat_black, 
  signaling = pathway_name, 
  color.heatmap = "Reds",
  title.name = paste(pathway_name, "signaling in Black"))

p_white <- netVisual_heatmap(
  cellchat_white, 
  signaling = pathway_name, 
  color.heatmap = "Reds",
  title.name = paste(pathway_name, "signaling in White"))

# --- Combine and Save the Plot ---
# (Remember: These are ComplexHeatmap objects, so we can't use ggsave)

combined_spp1_plot <- p_black + p_white


# --- Save to a file ---
png(filename = "f7d.png", 
    width = 12, 
    height = 6.1, 
    units = "in", 
    res = 300)

draw(combined_spp1_plot)

dev.off()
