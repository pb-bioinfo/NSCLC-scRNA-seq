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

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')



cellchat <- qread("cellchat_total.qs")
cellchat_black <- qread("cellchat_black.qs")
cellchat_white <- qread("cellchat_white.qs")

# the seurat object that contains only tumor cells
rds2 <- qread("rds2.qs")

#####

s7a <- DimPlot(rds2, reduction = "umap.rpca", label = F,  cols = col_vector, group.by = "stage2") +
  labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle("Functional Subtypes")

ggsave(filename = "s7a.png", plot = s7a, width = 12, height = 3.5, units = "in")





#####

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
s7b <- gg1 + gg2

#ComplexHeatmap, cannot use ggsave

png(filename = "s7b.png", 
    width = 12, 
    height = 6.1, 
    units = "in", 
    res = 300) 
draw(s7b)

dev.off()


#####

object.list <- list(B = cellchat_black, W = cellchat_white)

gg <- list()

manual_x_lim <- c(0, 120)

manual_y_lim <- c(0, 60)


for (i in 1:length(object.list)) {
  
  # Create the base plot (calculates math for the WHOLE network)
  net_plot <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    slot.name = "netP",
    title = names(object.list)[i],
    weight.MinMax = c(1000, 8000)
  )
  
  
  gg[[i]] <- net_plot +
    coord_cartesian(xlim = manual_x_lim, ylim = manual_y_lim) +
    theme(
      plot.title = element_text(size = 14)
    )+
    theme(axis.text.x = element_text( size =12),axis.title.x = element_text( size =12),
          axis.text.y = element_text(size = 12),axis.title.y = element_text( size =12))
  
}

s7c <- gg[[1]] | gg[[2]]

ggsave(filename = "s7c.png", plot = s7c, width = 12, height = 4.5, units = "in")




#####

object.list <- list(B = cellchat_black, W = cellchat_white)


total_types <- levels(cellchat@idents[[2]])

target_types <- total_types[c(21:25,29:41)]


gg <- list()

manual_x_lim <- c(0, 40)

manual_y_lim <- c(0, 60)


for (i in 1:length(object.list)) {
  # Create the base plot (calculates math for the WHOLE network)
  net_plot <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    slot.name = "netP",
    title = names(object.list)[i],
    weight.MinMax = c(1000, 6000)
  )
  
  
  net_plot$data <- net_plot$data[net_plot$data$labels %in% target_types, ]
  
  gg[[i]] <- net_plot +
    coord_cartesian(xlim = manual_x_lim, ylim = manual_y_lim) +
    theme(
      plot.title = element_text(size = 14)
    )+
    theme(axis.text.x = element_text( size =12),axis.title.x = element_text( size =12),
          axis.text.y = element_text(size = 12),axis.title.y = element_text( size =12))
  
}

s7d <- gg[[1]] | gg[[2]]


ggsave(filename = "s7d.png", plot = s7d, width = 12, height = 4.5, units = "in")


#####

#total incoming and outgoing interaction numbers 
extract_group_crosstalk_total <- function(obj, group1, group2) {
  count_mat <- obj@net$count
  weight_mat <- obj@net$weight
  
  
  g1_out <- rowSums(count_mat[group1, c(1:41), drop = FALSE])
  g1_in  <- colSums(count_mat[c(1:41), group1, drop = FALSE])
  g1_wt  <- rowSums(weight_mat[group1,  c(1:41), drop = FALSE]) + 
    colSums(weight_mat[ c(1:41), group1, drop = FALSE])
  
  df1 <- data.frame(labels = group1, Out_Count = g1_out, In_Count = g1_in, Total_Strength = g1_wt)
  
  g2_out <- rowSums(count_mat[group2,  c(1:41), drop = FALSE])
  g2_in  <- colSums(count_mat[ c(1:41), group2, drop = FALSE])
  g2_wt  <- rowSums(weight_mat[group2,  c(1:41), drop = FALSE]) + 
    colSums(weight_mat[ c(1:41), group2, drop = FALSE])
  
  df2 <- data.frame(labels = group2, Out_Count = g2_out, In_Count = g2_in, Total_Strength = g2_wt)
  
  return(rbind(df1, df2))
}



mye_cells <- levels(cellchat@idents[[1]])[21:25] 
t_cells <- levels(cellchat@idents[[1]])[29:39] 

df_black_counts_mye_T <- extract_group_crosstalk_total(cellchat_black, mye_cells, t_cells)
df_white_counts_mye_T <- extract_group_crosstalk_total(cellchat_white, mye_cells, t_cells)

df_fc <- merge(df_black_counts_mye_T, df_white_counts_mye_T, by = "labels", suffixes = c("_Black", "_White"))

df_fc <- df_fc %>%
  mutate(
    # Log2( (Black + 1) / (White + 1) ), not needed for the counts
    Log2FC_Out = log2((Out_Count_Black ) / (Out_Count_White )),
    Log2FC_In  = log2((In_Count_Black)  / (In_Count_White )),
    
    # Keep average strength for bubble size
    Avg_Strength = (Total_Strength_Black + Total_Strength_White) / 2
  )


df_fc_target <-df_fc

s7e <- ggplot(df_fc_target, aes(x = Log2FC_Out, y = Log2FC_In)) +
  # Crosshairs at 0 (0 means exactly equal between Black and White)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  
  geom_point(aes(size = Avg_Strength, color = labels), alpha = 0.8) +
  geom_text_repel(aes(label = labels), size = 5, box.padding = 0.5) +
  
  theme_bw()+labs(
    title = "Ratio Interaction Counts: Myeloid & T",
    x = "Log2 FC of Outgoing Counts (Black/White)",
    y = "Log2 FC of Incoming Counts (Black/White)",
    size = "Average Total\nStrength",
    color = "Cell Type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 12)
  )


ggsave(filename = "s7e.png", plot = s7e, width = 9, height = 9, units = "in")





##### 

extract_group_crosstalk <- function(obj, group1, group2) {
  count_mat <- obj@net$count
  weight_mat <- obj@net$weight
  
  
  g1_out <- rowSums(count_mat[group1, group2, drop = FALSE])
  g1_in  <- colSums(count_mat[group2, group1, drop = FALSE])
  g1_wt  <- rowSums(weight_mat[group1, group2, drop = FALSE]) + 
    colSums(weight_mat[group2, group1, drop = FALSE])
  
  df1 <- data.frame(labels = group1, Out_Count = g1_out, In_Count = g1_in, Total_Strength = g1_wt)
  
  g2_out <- rowSums(count_mat[group2, group1, drop = FALSE])
  g2_in  <- colSums(count_mat[group1, group2, drop = FALSE])
  g2_wt  <- rowSums(weight_mat[group2, group1, drop = FALSE]) + 
    colSums(weight_mat[group1, group2, drop = FALSE])
  
  df2 <- data.frame(labels = group2, Out_Count = g2_out, In_Count = g2_in, Total_Strength = g2_wt)
  
  return(rbind(df1, df2))
}



b_cells <- levels(cellchat@idents[[1]])[1:4] 
t_cells <- levels(cellchat@idents[[1]])[c(29:35,37:39)] 


df_black_targeted <- extract_group_crosstalk(cellchat_black, b_cells, t_cells)
df_white_targeted <- extract_group_crosstalk(cellchat_white, b_cells, t_cells)

df_diff_targeted <- merge(df_black_targeted, df_white_targeted, by = "labels", suffixes = c("_Black", "_White"))

df_diff_targeted <- df_diff_targeted %>%
  mutate(
    Delta_Out_Count = Out_Count_Black - Out_Count_White,
    Delta_In_Count  = In_Count_Black - In_Count_White,
    Avg_Strength    = (Total_Strength_Black + Total_Strength_White) / 2
  )

df_diff_targeted <- df_diff_targeted %>%
  mutate(
    Ratio_Out_Count = log2(Out_Count_Black / Out_Count_White),
    Ratio_In_Count  = log2(In_Count_Black / In_Count_White),
  )



s7f <- ggplot(df_diff_targeted, aes(x = Ratio_Out_Count, y = Ratio_In_Count)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  
  geom_point(aes(size = Avg_Strength, color = labels), alpha = 0.8) +  
  theme_bw() +labs(
    title = "Ratio Interaction Counts: B & T Only",
    x = "Log2 FC of Outgoing Counts (Black/White)",
    y = "Log2 FC of Incoming Counts (Black/White)",
    size = "Avg Interaction\nStrength",
    color = "Cell Type"
  )  +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 12)
  )

ggsave(filename = "s7f.png", plot = s7f, width = 9, height = 9, units = "in")





#####

mali_cells <- c(levels(cellchat@idents[[1]])[14])
cd8_cells <- levels(cellchat@idents[[1]])[36:39] 

nk_cells <- levels(cellchat@idents[[1]])[26:28] 
cyto_cells <- c(cd8_cells,nk_cells)




df_black_targeted_cyto_mali <- extract_group_crosstalk(cellchat_black, cyto_cells, mali_cells)
df_white_targeted_cyto_mali <- extract_group_crosstalk(cellchat_white, cyto_cells, mali_cells)

df_diff_targeted_cyto_mali <- merge(df_black_targeted_cyto_mali, df_white_targeted_cyto_mali, by = "labels", suffixes = c("_Black", "_White"))

df_diff_targeted_cyto_mali <- df_diff_targeted_cyto_mali %>%
  mutate(
    Delta_Out_Count = Out_Count_Black - Out_Count_White,
    Delta_In_Count  = In_Count_Black - In_Count_White,
    Avg_Strength    = (Total_Strength_Black + Total_Strength_White) / 2
  )

df_diff_targeted_cyto_mali <- df_diff_targeted_cyto_mali %>%
  mutate(
    Ratio_Out_Count = log2(Out_Count_Black / Out_Count_White),
    Ratio_In_Count  = log2(In_Count_Black / In_Count_White),
  )



s7g <- ggplot(df_diff_targeted_cyto_mali, aes(x = Ratio_Out_Count, y = Ratio_In_Count)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  
  geom_point(aes(size = Avg_Strength, color = labels), alpha = 0.8) +  
  theme_bw() +labs(
    title = "Ratio Interaction Counts: B & T Only",
    x = "Log2 FC of Outgoing Counts (Black/White)",
    y = "Log2 FC of Incoming Counts (Black/White)",
    size = "Avg Interaction\nStrength",
    color = "Cell Type"
  )  +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 12)
  )

ggsave(filename = "s7g.png", plot = s7g, width = 7, height = 9, units = "in")





#####


#Mono-macs \u2192 CD4 T Cells

cd4_t_cells <- c(29:35)

sources.use = 23
targets.use = cd4_t_cells



plot_max2 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 2,
                              angle.x = 45,
                              remove.isolate = TRUE)+
  labs(title = "Higher in White")


plot_max1 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 1,
                              angle.x = 45,
                              remove.isolate = TRUE)+
  labs(title = "Higher in Black")

s7h <- plot_max1 + plot_max2 +
  plot_annotation(
    title = "Mono-macs \u2192 CD4 T Cells",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
  )

ggsave(filename = "s7h.png", plot = s7h, width = 12, height = 9, units = "in")



#####


cd8_t_cells <- c(36:39)
sources.use = 24
targets.use = cd8_t_cells


plot_max2 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 2,
                              angle.x = 45,
                              remove.isolate = TRUE)+
  labs(title = "Higher in White")


plot_max1 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 1,
                              angle.x = 45,
                              remove.isolate = TRUE)+
  labs(title = "Higher in Black")

s7i <- plot_max1 + plot_max2 +
  plot_annotation(
    title = "Monocytes \u2192 CD8 T Cells",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
  )

ggsave(filename = "s7i.png", plot = s7i, width = 9, height = 7, units = "in")


#####



cd8_t_cells <- c(36:39)
#Mono-macs \u2192 CD8 T Cells
sources.use = 23
targets.use = cd8_t_cells


plot_max2 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 2,
                              angle.x = 45,
                              remove.isolate = TRUE)+
  labs(title = "Higher in White")


plot_max1 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 1,
                              angle.x = 45,
                              remove.isolate = TRUE)+
  labs(title = "Higher in Black")

s7j <- plot_max1 + plot_max2 +
  plot_annotation(
    title = "Mono-macs \u2192 CD8 T Cells",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
  )

ggsave(filename = "s7j.png", plot = s7j, width = 9, height = 7, units = "in")

#####




cd8_t_cells <- c(36:39)
#Mono-macs \u2192 CD8 T Cells
sources.use = 25
targets.use = cd8_t_cells


plot_max2 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 2,
                              angle.x = 45,
                              remove.isolate = TRUE)+
  labs(title = "Higher in White")


plot_max1 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 1,
                              angle.x = 45,
                              remove.isolate = TRUE)+
  labs(title = "Higher in Black")

s7k <- plot_max1 + plot_max2 +
  plot_annotation(
    title = "TAMs \u2192 CD8 T Cells",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
  )

ggsave(filename = "s7k.png", plot = s7k, width = 9, height = 7, units = "in")



#####

obj_T <- qread("obj_T.qs")

# 1 = Activated B
# 2 = Memory B
# 5 = Plasmablast
target_groups <- c(1,2,4)
names(target_groups) <- c("Activated B", "Memory B","Plasma")

# CD4 group: indices 29 to 35
cd4_t_cells <- c(29:35,40)
# CD8/Other group: indices 36 to 41
cd8_other_t_cells <- c(36,37,38,39,41)


group_name = names(target_groups)[1]

plot1 <- netVisual_bubble(cellchat,
                          sources.use = target_groups[1],
                          targets.use = cd4_t_cells,
                          comparison = c(1, 2),
                          angle.x = 45,
                          remove.isolate = TRUE)+ggtitle( paste0(group_name, " Cells \u2192 CD4+ T Cells"))


plot2 <- netVisual_bubble(cellchat,
                          sources.use = cd4_t_cells,
                          targets.use = target_groups[1],
                          comparison = c(1, 2),
                          angle.x = 45,
                          remove.isolate = TRUE)+ggtitle(paste0("CD4+ T Cells \u2192 ",group_name, " Cells"))

s7l <- plot1 | plot2

ggsave(filename = "s7l.png", plot = s7l, width = 12, height = 8, units = "in")



#####




group_name = names(target_groups)[2]

plot1 <- netVisual_bubble(cellchat,
                          sources.use = target_groups[2],
                          targets.use = cd8_other_t_cells,
                          comparison = c(1, 2),
                          angle.x = 45,
                          remove.isolate = TRUE)+ggtitle( paste0(group_name, " Cells \u2192 Effector T Cells"))


plot2 <- netVisual_bubble(cellchat,
                          sources.use = cd8_other_t_cells,
                          targets.use = target_groups[2],
                          comparison = c(1, 2),
                          angle.x = 45,
                          remove.isolate = TRUE)+ggtitle(paste0("Effector T Cells \u2192 ",group_name, " Cells"))

s7m <- plot1 | plot2

ggsave(filename = "s7m.png", plot = s7m, width = 12, height = 8, units = "in")



#####

group_name = names(target_groups)[3]

plot1 <- netVisual_bubble(cellchat,
                          sources.use = target_groups[3],
                          targets.use = cd8_other_t_cells,
                          comparison = c(1, 2),
                          angle.x = 45,
                          remove.isolate = TRUE)+ggtitle( paste0(group_name, " Cells \u2192 Effector T Cells"))


plot2 <- netVisual_bubble(cellchat,
                          sources.use = cd8_other_t_cells,
                          targets.use = target_groups[3],
                          comparison = c(1, 2),
                          angle.x = 45,
                          remove.isolate = TRUE)+ggtitle(paste0("Effector T Cells \u2192 ",group_name, " Cells"))

s7n <- plot1 | plot2

ggsave(filename = "s7n.png", plot = s7n, width = 12, height = 8, units = "in")



#####

tls_sigs <- c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10","CXCL11","CXCL13","CXCL3")


sce_score <- AddModuleScore(
  object = obj_T,
  features = list(tls_sigs),
  name = 'score_TLS'
)





meta_data <- as.data.frame(sce_score@meta.data)

s7o <- ggplot(meta_data, aes(x = score_TLS1, colour = R_T_N)) +scale_color_manual(values = color.liberal)+
  stat_ecdf(geom = "step", size = 0.5) + # Use stat_ecdf with geom = "step"
  labs(
    title = "",
    x = "TLS Score",
    y = "Cumulative Proportion",
    colour = " "
  )+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size =10),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 10))+
  theme(axis.line = element_line(colour = "black"))


ggsave(filename = "s7o.png", plot = s7o, width = 4.5, height = 3, units = "in",dpi = 300)



#####

cd8_other_t_cells <- c(36,37,38,39,41)
nk_cells <- c(26:28) 



sources.use = cd8_other_t_cells
targets.use = 14


plot_max2 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 2,
                              angle.x = 45,
                              remove.isolate = FALSE)+
  labs(title = "Higher in White")


plot_max1 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 1,
                              angle.x = 45,
                              remove.isolate = FALSE)+
  labs(title = "Higher in Black")

s7p <- plot_max1 + plot_max2 +
  plot_annotation(
    title = "Cytotoxic/Effector T Cells \u2192 Malignant Epithelial Cells",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
  )

ggsave(filename = "s7p.png", plot = s7p, width = 10, height = 9, units = "in")


#####


sources.use = nk_cells
targets.use = 14


plot_max2 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 2,
                              angle.x = 45,
                              remove.isolate = FALSE)+
  labs(title = "Higher in White")


plot_max1 <- netVisual_bubble(cellchat,
                              sources.use = sources.use,
                              targets.use = targets.use,
                              comparison = c(1, 2),
                              max.dataset = 1,
                              angle.x = 45,
                              remove.isolate = FALSE)+
  labs(title = "Higher in Black")

s7q <- plot_max1 + plot_max2 +
  plot_annotation(
    title = "NK Cells \u2192 Malignant Epithelial Cells",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
  )

ggsave(filename = "s7q.png", plot = s7q, width = 10, height = 9, units = "in")

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
png(filename = "s7r.png", 
    width = 12, 
    height = 6.1, 
    units = "in", 
    res = 300)

draw(combined_spp1_plot)

dev.off()

