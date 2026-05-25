library(scales)  # For DiscretePalette function
library(grid)
library(ggplot2)
library(ggtext)
library(purrr)
library(qs)
library(ggpubr)
library(hrbrthemes)
library(patchwork)
library(cowplot)   # For get_legend() and ggdraw()

library(data.table)
library(reshape2)
library(viridis)
library(ggsci)
library(qusage)
library(tidyverse)
library(RColorBrewer)

library(DescTools)

library(stats)
library(lsa)
library(dplyr)


library(Seurat)
library(igraph)
library(presto)


setwd("~/diskD/lung/paneltotal/fig4/")


obj_myeloid <-qread("obj_myeloid.qs")
obj_macro <-qread("obj_macro.qs")
obj_DC <-qread("obj_DC.qs")

#the subset of cells that monocle2 analysis performed on
trajectory_subset <- qread("trajectory_subset.qs")



#monocyte
obj_mono <-qread("obj_mono.qs")

#mono-mac
obj_monomac <-qread("obj_monomac.qs")



color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#####




s4a <- DimPlot(obj_myeloid, reduction = "umap.Myeloid.rpca", label=T, cols = col_vector,group.by = "seurat_clusters",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")
s4a

ggsave(filename = "s4a.tiff", plot = s4a, width = 4.5, height = 4, units = "in")


#####


smalllist = c("AIF1","FCN1","MAF","APOE","FABP4","LAMP3","CCL17")

plots <- lapply(1:7, function(i) {
  if (i == 4){
    FeaturePlot(obj_myeloid,features = smalllist[i], reduction = "umap.Myeloid.rpca",cols = c('grey','red'))+labs(y = "UMAP_2", x = "UMAP_1",max.cutoff = 4.6)+theme(legend.position = "none")
  }else{
    FeaturePlot(obj_myeloid,features = smalllist[i], reduction = "umap.Myeloid.rpca",cols = c('grey','red'))+labs(y = " ", x = " ",max.cutoff = 4.6,no.axes = TRUE)+theme(legend.position = "none")
  }
})
s4b <- wrap_plots(plots, ncol = 4) +  plot_layout(guides = 'collect') +
  guides(fill = guide_legend(title = "Expression")) +
  theme(legend.position = "right")  

ggsave(filename = "s4b.tiff", plot = s4b, width = 7, height = 3.6,units = "in")

#####


genes_self_TAM <- c("SPP1", "RNASE1", "GPNMB","SELENOP","VEGFA", "SLAMF9")
genes_self_cDC <- c("FCER1A", "CCL17", "IDO1", "CXCL9", "CLEC4C", "CCR7", "LAMP3")
genes_self_ave <- c("FABP4", "CES1", "MARCO", "MCEMP1", "PPARG")
genes_self_prolif <- unlist(strsplit("STMN1, H2AFZ, TUBA1B, PCNA", ", "))
genes_self_proinflammatory <- c("CCL15","CXCL13")
genes_NatCom2020_monocypte <- unlist(strsplit("CTSS, FCN1, S100A8, S100A9, LYZ, VCAN, PLAC8", ", "))
genes_NatCom2020_anti_inflammatory <- unlist(strsplit("APOE, SELENOP, C1QA, C1QB, C1QC", ", "))
genes_NatCom2020_mo_mac <- unlist(strsplit("MAFB, MAF, CX3CR1, ITGAM, CSF1R", ", "))
genes_NatCom2020_macrophage <- unlist(strsplit("LGMN, CTSB, CD14, FCGR3A", ", "))
genes_monoDC <- c("CD14","FCGR1B", "CLEC10A", "MRC1")
genes_pDC <- c("IRF4","LILRA4","TCF4")
genes_granulocyte <- c("G0S2", "S100A12", "FCGR3B", "IFITM3", "FCER1A", "S100A8", "S100A9")

my_levels <- paste0("g", 
                    c(4,7,17,13,
                      11,14,3,29,
                      10,8,18,9,5,27,32,1,21,
                      31,26,0,15,6,19,22,28,16,30,12,
                      20,25,33,2,23,24
                    ))

markers3 <- c(genes_NatCom2020_monocypte,
              genes_NatCom2020_macrophage,
              genes_self_TAM,
              genes_NatCom2020_mo_mac,
              genes_self_cDC,genes_pDC,genes_self_ave,genes_self_prolif,genes_self_proinflammatory,genes_NatCom2020_anti_inflammatory,genes_granulocyte)


Idents(obj_Myeloid) <- "RNA_snn_res.1"

cluster.averages = AverageExpression(obj_Myeloid, return.seurat = T, group.by = 'RNA_snn_res.1')

Idents(cluster.averages) <- factor(cluster.averages$seurat_clusters, levels = my_levels)

s4c <- DoHeatmap(
  object = cluster.averages,
  features = markers3,  # Genes
  size = 5,
  draw.lines = FALSE,
  disp.min = 0,
  disp.max = 2,
) +
  scale_fill_gradientn(colors = c("lightgrey", "blue")) +
  #scale_x_discrete(labels = gsub("^g", "", levels(cluster.averages$seurat_clusters_clean))) +  
  theme(
    axis.text.y = element_text(size = 10, family = "Arial", color = "black"),  
    axis.text.x.bottom = element_text(size = 8, angle = 90,hjust = 1, family = "Arial", color = "black"), 
    axis.text.x.top = element_blank(),  #axis.text.x.top = element_text(size = 8, angle = 45,  family = "Arial",hjust = 1, family = "Arial", color = "black"),  # Top axis cluster labels
    legend.position = "none", 
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank()   # Remove y-axis title
  )


ggsave(filename = "s4c.tiff", plot = s4c, width = 6.5, height = 9, units = "in")


#####


s4d <- DimPlot(obj_myeloid, reduction = "umap.Myeloid.rpca", label=T, cols = color.liberal,group.by = "stage1_v3",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")
ggsave(filename = "s4d.tiff", plot = s4d, width = 3.4, height = 3, units = "in")


#####


pt <- table(obj_myeloid@meta.data$Sample_ID, obj_myeloid@meta.data$stage1_v3)
pt <- as.data.frame(pt)




s4e1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x =element_blank())+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s4e2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))







combined_plot <- (s4e1 / s4e2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot

ggsave(filename = "s4e.tiff", plot = combined_plot, width = 11, height = 5, units = "in")

#####


s4f <- DimPlot(obj_myeloid, reduction = "umap.Myeloid.rpca",split.by = "T_N" , label=F,group.by ="stage1_v3",ncol = 2)+ NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = color.liberal)+
  theme(strip.text = element_text(face = "bold", size = 15))

ggsave(filename = "s4f.tiff", plot = s4f, width = 5.6, height = 3, units = "in")


#####

pt <- table(obj_myeloid@meta.data$T_N, obj_myeloid@meta.data$stage1_v3)
pt <- as.data.frame(pt)

s4g1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))+NoLegend()


s4g2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))


combined_plot <- (s4g1 | s4g2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot

ggsave(filename = "s4g.tiff", plot = combined_plot, width = 4.3, height = 3, units = "in")


#####



custom_colors <- c("orange",      # orange
                   "#377EB8",      # blue
                   "#A50F15",      # dark red
                   "#9ECAE1")



# -----------------------------
# 1. Build patient-level proportions
# -----------------------------
meta_df <- obj_T@meta.data %>%
  dplyr::select(
    obj_myeloid,
    celltype = "stage1_v3",
    group = "T_N"
  ) 

celltype_levels <- unique(meta_df$celltype)
total_rows <- floor(length(celltype_levels) / 3)

sample_group_df <- meta_df %>%
  dplyr::distinct(Sample_ID, group)

group_df_long <- meta_df %>%
  dplyr::count(Sample_ID, group, celltype, name = "n") %>%
  dplyr::right_join(
    tidyr::expand_grid(
      sample_group_df,
      celltype = celltype_levels
    ),
    by = c("Sample_ID", "group", "celltype")
  ) %>%
  dplyr::mutate(n = tidyr::replace_na(n, 0)) %>%
  dplyr::group_by(Sample_ID, group) %>%
  dplyr::mutate(value = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    Var3 = group,
    variable = celltype,
    value = value
  )

total_types <- as.character(unique(group_df_long$variable))

# -----------------------------
# 2. One-sided Wilcoxon tests
# -----------------------------
comparisons_less <- ggpubr::compare_means(
  value ~ Var3,
  group.by = "variable",
  data = group_df_long,
  method = "wilcox.test",
  alternative = "less"
)

comparisons_greater <- ggpubr::compare_means(
  value ~ Var3,
  group.by = "variable",
  data = group_df_long,
  method = "wilcox.test",
  alternative = "greater"
)

combined_significant_comparisons <- dplyr::bind_rows(
  comparisons_greater %>%
    dplyr::filter(p < 0.05) %>%
    dplyr::mutate(alternative = "greater"),
  
  comparisons_less %>%
    dplyr::filter(p < 0.05) %>%
    dplyr::mutate(alternative = "less")
)


# -----------------------------
# 3. Plot one cell type
# -----------------------------
plot_one_type <- function(i) {
  
  subdata <- group_df_long %>%
    dplyr::filter(variable == total_types[i]) 
  
  maxx <- max(subdata$value, na.rm = TRUE)
  
  onetype_plot <- ggplot(
    subdata,
    aes(x = Var3, y = value, fill = Var3)
  ) +
    geom_boxplot() +
    hrbrthemes::theme_ipsum() +
    facet_wrap(
      ~ variable,
      scales = "free",
      labeller = label_wrap_gen(width = 20)
    ) +
    scale_fill_manual(values = custom_colors) +
    ylim(0, maxx + 0.4) +
    geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    labs(x = "", y = "Proportion", fill = "") +
    theme(
      strip.text = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")
    ) +
    Seurat::NoLegend()
  
  if (i == 1) {
    onetype_plot <- onetype_plot +
      labs(y = "Proportion") +
      theme(axis.title.y = element_text(size = 15))
  } else {
    onetype_plot <- onetype_plot +
      labs(y = "")
  }
  
  filtered_comparisons <- combined_significant_comparisons %>%
    dplyr::filter(variable == total_types[i])
  
  if (nrow(filtered_comparisons) > 0) {
    for (g in seq_len(nrow(filtered_comparisons))) {
      
      pair <- c(filtered_comparisons$group1[g], filtered_comparisons$group2[g])
      
      
      plot_pair <- list(rev(pair))
      
      onetype_plot <- onetype_plot +
        ggpubr::stat_compare_means(
          method = "wilcox.test",
          comparisons = plot_pair,
          size = 6,
          label = "p.signif",
          hide.ns = TRUE,
          data = subdata,
          method.args = list(alternative = filtered_comparisons$alternative[g]),
          label.y = maxx + g * 0.05
        )
      
    }
  }
  
  onetype_plot <- onetype_plot +
    theme(axis.text.x = element_blank())
  
  return(onetype_plot)
}

# -----------------------------
# 4. Generate plots and shared legend
# -----------------------------
plot_list <- lapply(seq_along(total_types), plot_one_type)

legend_plot <- plot_list[[1]] +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 4))

shared_legend_grob <- cowplot::get_legend(legend_plot)
legend_for_grid <- cowplot::ggdraw(shared_legend_grob)

plot_list_with_legend <- c(plot_list, list(legend_for_grid))

final_plot <- patchwork::wrap_plots(
  plot_list_with_legend,
  ncol = 3
) +
  patchwork::plot_layout(guides = "keep") +
  patchwork::plot_annotation() &
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5),
    panel.spacing = unit(0.1, "lines")
  )


s4h <- final_plot

ggsave(filename = "s4h.tiff", plot = s4h, width = 4.5, height = 4, units = "in")


#####


Idents(obj_macro)<- "Race_TNM"

s4i <- DEenrichRPlot(
  obj_macro,
  ident.1 = "B_T",
  ident.2 = "W_T",
  balanced = TRUE,
  logfc.threshold = 0.25, # Lowered slightly to catch more genes for pathway analysis
  assay = "RNA",         # Usually safer to specify "RNA"
  max.genes = 100,       # Increased: KEGG needs more genes than Hallmark to trigger hits
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "MSigDB_Hallmark_2020",
  num.pathway = 10,      # Increased to see more results
  return.gene.list = FALSE)


s4i[[1]] <- s4i[[1]]  + 
  labs(title =paste0("Black Tumor ","TAM"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


s4i[[2]] <- s4i[[2]]  + 
  labs(title =paste0("White Tumor ","TAM"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


ggsave(filename = "s4i.tiff", plot = s4i, width = 6.2, height = 3, units = "in")

#####
antigen_presentation_genes <- c(
  "CD1B", "CD1D", "CD40", "CD80", "CD86",
  "HLA-A", "HLA-B", "HLA-C",
  "TAP1", "TAP2", "ERAP1", "ERAP2",   
  "HLA-DRA", "HLA-DRB1", "HLA-DMA", "HLA-DMB", "HLA-DQA1", "HLA-DQB1"
)


Idents(obj_macro) <- "Type2"

markers <- FindMarkers(
  obj_macro,
  ident.1 = "LUAD",
  ident.2 = "LUSC",
  logfc.threshold = logfc_threshold,
  min.pct = min_pct,
  features = antigen_presentation_genes 
)


subset_degs <- markers[rownames(markers) %in% antigen_presentation_genes, ]

if (nrow(subset_degs) == 0) {
  warning("No significant markers found from the provided gene list.")
}

subset_degs$gene <- rownames(subset_degs)

print(subset_degs)

# 3. Create Plot
s4j <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw() +
  labs(
    title = "Enriched in LUAD TAMs",
    x = "",
    y = "Average Log2FC",
    size = "-log10(p-adj)",
    color = "Log2FC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.box = "horizontal")+
  guides( color = "none")



ggsave(filename = "s4j.tiff", plot = s4j, width = 5, height = 4, units = "in")

#####



Idents(obj_monomac) <- "Type2"

markers <- FindMarkers(
  obj_monomac,
  ident.1 = "LUAD",
  ident.2 = "LUSC",
  logfc.threshold = logfc_threshold,
  min.pct = min_pct,
  features = antigen_presentation_genes 
)


subset_degs <- markers[rownames(markers) %in% antigen_presentation_genes, ]

if (nrow(subset_degs) == 0) {
  warning("No significant markers found from the provided gene list.")
}

subset_degs$gene <- rownames(subset_degs)

print(subset_degs)

# 3. Create Plot
s4j2 <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw() +
  labs(
    title = "Enriched in LUAD Mono-Macs",
    x = "",
    y = "Average Log2FC",
    size = "-log10(p-adj)",
    color = "Log2FC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.box = "horizontal")+
  guides( color = "none")



ggsave(filename = "s4j2.tiff", plot = s4j2, width = 5, height = 4, units = "in")


#####


Idents(obj_mono) <-  "Race_TNM"

# Here we run standard DE to get adjusted p-values properly.
markers <- FindMarkers(
  obj_mono,
  ident.1 = "B_T",
  ident.2 = "W_T",
  logfc.threshold = 0.1,
  min.pct = 0.1,
  features = antigen_presentation_genes 
)

# Filter to ensure we only have our genes
subset_degs <- markers[rownames(markers) %in% antigen_presentation_genes, ]

if (nrow(subset_degs) == 0) {
  warning("No significant markers found from the provided gene list.")
}

subset_degs$gene <- rownames(subset_degs)

# Print the table for verification
print(subset_degs)

# 3. Create Plot
s4l <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw() +
  labs(
    title = "Enriched in black tumor Monocytes",
    x = "",
    y = "Average Log2FC",
    size = "-log10(p-adj)",
    color = "Log2FC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.box = "horizontal")+
  guides( color = "none")


ggsave(filename = "s4l.tiff", plot = s4l, width = 4.5, height = 3.5, units = "in")


#####

Idents(obj_monomac) <-  "Race_TNM"

# Here we run standard DE to get adjusted p-values properly.
markers <- FindMarkers(
  obj_monomac,
  ident.1 = "B_T",
  ident.2 = "W_T",
  logfc.threshold = 0.1,
  min.pct = 0.1,
  features = antigen_presentation_genes 
)

# Filter to ensure we only have our genes
subset_degs <- markers[rownames(markers) %in% antigen_presentation_genes, ]

if (nrow(subset_degs) == 0) {
  warning("No significant markers found from the provided gene list.")
}

subset_degs$gene <- rownames(subset_degs)

# Print the table for verification
print(subset_degs)

# 3. Create Plot
s4l2 <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw() +
  labs(
    title = "Enriched in black tumor Mono-Macs",
    x = "",
    y = "Average Log2FC",
    size = "-log10(p-adj)",
    color = "Log2FC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.box = "horizontal")+
  guides( color = "none")

ggsave(filename = "s4l2.tiff", plot = s4l2, width = 4.5, height = 3.5, units = "in")


#####

traj_mono = subset(trajectory_subset,stage1_res.1_v2 == "Monocyte")

Idents(traj_mono) <- "State"

s4m <- DEenrichRPlot( 
  traj_mono,
  ident.1 = 9,
  ident.2 = 6,
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = NULL,
  max.genes=100,
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = "MSigDB_Hallmark_2020",
  num.pathway = 10,
  return.gene.list = FALSE)


s4m[[1]] <- s4m[[1]]  + 
  labs(title =paste0("State 9 ","Monocyte"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


s4m[[2]] <- s4m[[2]]  + 
  labs(title =paste0("State 6 ","Monocyte"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))




ggsave(filename = "s4m.tiff", plot = s4m, width = 6.2, height = 3, units = "in")



#####

traj_macro = subset(trajectory_subset,stage1_v3 == "TAM")

Idents(traj_macro) <- "State"

s4n <- DEenrichRPlot(
  traj_macro,
  ident.1 = 3,
  ident.2 = 5,
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = NULL,
  max.genes=100,
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = "MSigDB_Hallmark_2020",
  num.pathway = 10,
  return.gene.list = FALSE)

s4n[[1]] <- s4n[[1]]  + 
  labs(title =paste0("State 3 ","TAM"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


s4n[[2]] <- s4n[[2]]  + 
  labs(title =paste0("State 5 ","TAM"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))




ggsave(filename = "s4n.tiff", plot = s4n, width = 6.2, height = 3, units = "in")

#####


traj_monomac = subset(trajectory_subset,stage1_v3 == "Mono-mac")
Idents(traj_monomac) <- "State"

s4o <- DEenrichRPlot(
  traj_monomac,
  ident.1 = 3,
  ident.2 = 5,
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = NULL,
  max.genes=100,
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = "MSigDB_Hallmark_2020",
  num.pathway = 10,
  return.gene.list = FALSE)

s4o[[1]] <- s4o[[1]]  + 
  labs(title =paste0("State 3 ","Mono-mac"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


s4o[[2]] <- s4o[[2]]  + 
  labs(title =paste0("State 5 ","Mono-mac"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


ggsave(filename = "s4o.tiff", plot = s4o, width = 6.2, height = 3, units = "in")


#####

s4p <- DimPlot(obj_DC, reduction = "umap", label=T, cols = color.liberal,group.by = "seurat_clusters",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")

ggsave(filename = "s4p.tiff", plot = s4p, width = 4, height = 3.5, units = "in")

#####


s4q <- DimPlot(obj_DC, reduction = "umap", label=T, cols = color.liberal,group.by = "DC",label.size = 4,repel = TRUE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")

ggsave(filename = "s4q.tiff", plot = s4q, width = 4, height = 3.5, units = "in")



