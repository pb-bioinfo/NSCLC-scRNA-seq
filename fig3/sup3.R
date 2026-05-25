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


setwd("~/diskD/lung/paneltotal/fig3/")



color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f68231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')


obj_Fibro <-qread("obj_Fibro.qs")

obj_caf2 <- qread("obj_caf2.qs")

obj_Endo <- qread("obj_Endo.qs")


#####


s6a <- DimPlot(obj_Fibro, reduction = "umap.endo.rpca", label=T, cols = col_vector,group.by = "seurat_clusters",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")



ggsave(filename = "s6a.png", plot = s6a, width = 3.4, height = 3, units = "in")


#####

fibro_markers <- c(
  "RGCC", "MACF1", "A2M", "DST", "NPNT", #normal
  "ACTA2", "PTN", "MYL9", "FAM162B", "KCNK3", #MF
  "FAP", "COL3A1", "COL1A1", "COL1A2", "POSTN", #CAF
  "FGFBP1", "GPX2", "NPSR1", #CAF
  "CFD", "CXCL14", "OGN", "GSN" #immune-modulator
)

cluster.averages <- AverageExpression(obj_Fibro,group.by = "seurat_clusters",return.seurat = TRUE)

my_levels <-  c("g2", "g1", "g7", "g10","g0",   "g4", "g3", "g5", "g6", "g8","g9", "g11")

Idents(cluster.averages) <- factor(cluster.averages$seurat_clusters, levels = my_levels)


s3b <- DoHeatmap(
  object = cluster.averages,
  features = fibro_markers,  # Genes
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


ggsave(filename = "s3b.tiff", plot = s3b, width = 6.5, height = 9, units = "in")


#####

s3c <- DimPlot(obj_Fibro, reduction = "umap.endo.rpca", label=T, cols = color.liberal,group.by = "endo",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")

ggsave(filename = "s3c.png", plot = s3c, width = 3.4, height = 3, units = "in")



#####




pt <- table(obj_Fibro@meta.data$Sample_ID, obj_Fibro@meta.data$endo)
pt <- as.data.frame(pt)




s3d1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x =element_blank())+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s3d2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))







combined_plot <- (s3d1 / s3d2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")

ggsave(filename = "s3d.png", plot = combined_plot,  width = 11, height = 5, units = "in")


#####



pt <- table(obj_Fibro@meta.data$T_N, obj_Fibro@meta.data$endo)
pt <- as.data.frame(pt)




s3e1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s3e2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))






combined_plot <- (s3e1 | s3e2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")

ggsave(filename = "s3e.png", plot = combined_plot, width = 4.2, height = 3,, units = "in")



#####


obj_mf <-subset(obj_Fibro,endo == "MF")
Idents(obj_mf) <- "Race_TNM"

s3f <- DEenrichRPlot(
  obj_mf,
  ident.1 = "B_T",
  ident.2 = "W_T",
  balanced = TRUE,
  logfc.threshold = 0.1, # Lowered slightly to catch more genes for pathway analysis
  assay = "RNA",         # Usually safer to specify "RNA"
  max.genes = 500,       # Increased: KEGG needs more genes than Hallmark to trigger hits
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "KEGG_2021_Human", # <--- THIS IS THE CHANGE
  num.pathway = 10,      # Increased to see more results
  return.gene.list = FALSE)
s3f[[1]] <- s3f[[1]] + 
  labs(title =paste0("Black Tumor ","MF"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.75) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )


s3f[[2]] <- s3f[[2]] + 
  labs(title =paste0("White Tumor ","MF"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.75) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )


ggsave(filename = "s3f.png", plot = s3f, width = 9.6, height = 3, units = "in")

#####

obj_caf1 <- qread("obj_caf1.qs")

Idents(obj_caf1) <- "Race_TNM"

s3g <- DEenrichRPlot(
  obj_caf1,
  ident.1 = "B_T",
  ident.2 = "W_T",
  balanced = TRUE,
  logfc.threshold = 0.1, # Lowered slightly to catch more genes for pathway analysis
  assay = "RNA",         # Usually safer to specify "RNA"
  max.genes = 500,       # Increased: KEGG needs more genes than Hallmark to trigger hits
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "KEGG_2021_Human", # <--- THIS IS THE CHANGE
  num.pathway = 10,      # Increased to see more results
  return.gene.list = FALSE)


s3g[[1]] <- s3g[[1]] + 
  labs(title =paste0("Black Tumor ","CAF-1"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.75) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )


s3g[[2]] <- s3g[[2]] + 
  labs(title =paste0("White Tumor ","CAF-1"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.75) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )



ggsave(filename = "s3g.png", plot = s3g, width = 9.6, height = 3, units = "in")

#####



inflammatory_caf_genes <- c(
  "IL6",
  "IL11",
  "IL18",
  "LIF",
  "CSF2",  # This is the gene for GM-CSF
  "CXCL1",
  "CXCL2",
  "CXCL8",
  "CXCL10",
  "CCL2",
  "CCL8"
)



Idents(obj_caf2) <- "Race_TNM"

# Here we run standard DE to get adjusted p-values properly.
markers <- FindMarkers(
  obj_caf2,
  ident.1 = "B_T",
  ident.2 = "W_T",
  logfc.threshold = 0.1,
  min.pct = 0.1,
  features = inflammatory_caf_genes 
)

# Filter to ensure we only have our genes
subset_degs <- markers[rownames(markers) %in% inflammatory_caf_genes, ]

if (nrow(subset_degs) == 0) {
  warning("No significant markers found from the provided gene list.")
}

subset_degs$gene <- rownames(subset_degs)

# 3. Create Plot
s3h <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw() +
  labs(
    title = "Black Tumor vs White Tumor CAF-2",
    x = "",
    y = "Average Log2FC",
    size = "-log10(p-adj)",
    color = "Log2FC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.box = "horizontal")+
  guides( color = "none")

ggsave(filename = "s3h.png", plot = s3h, width = 4.3, height = 3.5, units = "in")


#####



s3i <- DimPlot(obj_Endo, reduction = "umap.endo.rpca", label=T, cols = col_vector,group.by = "seurat_clusters",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")

ggsave(filename = "s3i.png", plot = s3i, width = 3.4, height = 3, units = "in")


#####


markers3  <- c(
  "HSPG2", #Tumor EC
  "INSR",
  "VWA1", #Tumor EC
  "RGCC", #Tip-like
  "RAMP3",
  "ADM", #Tip-like
  "ACKR1", #Stalk
  "SELP",
  "GJA5", #Arterial
  "FBLN5",
  "TYROBP", #EPC
  "C1QB",
  "PDPN", #lymphatic
  "CCL21",
  "LYVE1",
  "PROX1",
  "FLT1", #VEGF Pathway
  "KDR",
  "FLT4",
  "VEGFA",
  "DLL1", #Notch Pathway
  "DLL4",
  "JAG1",
  "JAG2",
  "NOTCH1",
  "NOTCH4"
)
cluster.averages <- AverageExpression(obj_Endo,group.by = "seurat_clusters",return.seurat = TRUE)

my_levels <- as.vector(c("g3","g1","g2","g7","g10","g11","g13","g8","g12","g4","g9","g0","g6","g5"))

Idents(cluster.averages) <- factor(cluster.averages$seurat_clusters, levels = my_levels)


s3j <- DoHeatmap(
  object = cluster.averages,
  features = fibro_markers,  # Genes
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


ggsave(filename = "s3j.png", plot = s3j, width = 6.5, height = 9, units = "in")


#####
s3k <- DimPlot(obj_Endo, reduction = "umap.endo.rpca", label=T, cols = color.liberal,group.by = "endo",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")


ggsave(filename = "s3k.png", plot = s3k, width = 3.4, height = 3, units = "in")


#####



pt <- table(obj_Endo@meta.data$Sample_ID, obj_Endo@meta.data$endo)
pt <- as.data.frame(pt)




s3l1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x =element_blank())+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s3l2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))







combined_plot <- (s3l1 / s3l2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")

ggsave(filename = "s3l.png", plot = combined_plot,  width = 11, height = 5, units = "in")



#####


pt <- table(obj_Endo@meta.data$T_N, obj_Endo@meta.data$endo)
pt <- as.data.frame(pt)




s3m1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s3m2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))






combined_plot <- (s3m1 | s3m2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")

ggsave(filename = "s3m.png", plot = combined_plot, width = 4.2, height = 3,, units = "in")



#####




orders2 = c("Black Normal","White Normal", "Black Tumor","White Tumor")

custom_colors <- c(
  "Black Normal" = "orange",
  "White Normal" = "#377EB8",
  "Black Tumor"  = "#A50F15",
  "White Tumor"  = "#9ECAE1"
)

NOT_comparisons <- list(
  c("Black Normal", "White Tumor"),
  c("White Tumor", "Black Normal"),
  c("Black Tumor", "White Normal"),
  c("White Normal", "Black Tumor")
)

# -----------------------------
# 1. Build patient-level proportions
# -----------------------------
meta_df <- obj_Endo@meta.data %>%
  dplyr::select(
    Sample_ID,
    celltype = "endo",
    group = "R_T_N"
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

is_forbidden_comparison <- function(pair) {
  any(vapply(
    NOT_comparisons,
    function(x) identical(as.character(pair), as.character(x)),
    logical(1)
  ))
}

# -----------------------------
# 3. Plot one cell type
# -----------------------------
plot_one_type <- function(i) {
  
  subdata <- group_df_long %>%
    dplyr::filter(variable == total_types[i]) %>%
    dplyr::mutate(
      Group_Label = factor(Var3, levels = orders2)
    )
  
  maxx <- max(subdata$value, na.rm = TRUE)
  
  onetype_plot <- ggplot(
    subdata,
    aes(x = Group_Label, y = value, fill = Group_Label)
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
    geom_jitter(color = "black", size = 0.4, alpha = 0.9,height = 0) +
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
      
      if (!is_forbidden_comparison(pair)) {
        
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

s3n <- patchwork::wrap_plots(
  plot_list_with_legend,
  ncol = 3
) +
  patchwork::plot_layout(guides = "keep") +
  patchwork::plot_annotation() &
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5),
    panel.spacing = unit(0.1, "lines")
  )


ggsave(filename = "s3n.png", plot = s3n, width = 6, height = 7, units = "in")


#####


obj_TEC <- subset(obj_Endo,endo == "Tumor EC")

Idents(obj_TEC) <- "Race_TNM"

s3o <- DEenrichRPlot(
  obj_TEC,
  ident.1 = "B_T",
  ident.2 = "W_T",
  balanced = TRUE,
  logfc.threshold = 0.25,  
  assay = "RNA",          
  max.genes = 100,        
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "MSigDB_Hallmark_2020",  
  num.pathway = 10,      # Increased to see more results
  return.gene.list = FALSE)


s3o[[1]] <- s3o[[1]] + 
  labs(title =paste0("Black Tumor ","Tumor EC"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.9) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )


s3o[[2]] <- s3o[[2]] + 
  labs(title =paste0("White Tumor ","Tumor EC"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.9) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )
ggsave(filename = "s3o.png", plot = s3o, width = 7.3, height = 3, units = "in")


#####

obj_Pan_V <- subset(obj_Endo,!(endo %in%  c("Lymphatic","Unknown")))

Idents(obj_Pan_V) <- "Race_TNM"

# Here we run standard DE to get adjusted p-values properly.
markers <- FindMarkers(
  obj_Pan_V,
  ident.1 = "B_T",
  ident.2 = "W_T",
  logfc.threshold = 0.1,
  min.pct = 0.1,
  features = all_antigen_presentation_genes 
)

# Filter to ensure we only have our genes
subset_degs <- markers[rownames(markers) %in% all_antigen_presentation_genes, ]

if (nrow(subset_degs) == 0) {
  warning("No significant markers found from the provided gene list.")
}

subset_degs$gene <- rownames(subset_degs)

# 3. Create Plot
s3p <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw() +
  labs(
    title = "Black Tumor vs White Tumor Pan-Vascular EC",
    x = "",
    y = "Average Log2FC",
    size = "-log10(p-adj)",
    color = "Log2FC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.box = "horizontal")+
  guides( color = "none")

ggsave(filename = "s3p.png", plot = s3p, width = 4.3, height = 3.5, units = "in")


#####


orders2 = c("Black Normal","White Normal", "Black Tumor","White Tumor")

custom_colors <- c(
  "Black Normal" = "orange",
  "White Normal" = "#377EB8",
  "Black Tumor"  = "#A50F15",
  "White Tumor"  = "#9ECAE1"
)

NOT_comparisons <- list(
  c("Black Normal", "White Tumor"),
  c("White Tumor", "Black Normal"),
  c("Black Tumor", "White Normal"),
  c("White Normal", "Black Tumor")
)

# -----------------------------
# 1. Build patient-level proportions
# -----------------------------
meta_df <- obj_Endo@meta.data %>%
  dplyr::select(
    Sample_ID,
    celltype = "Angiogenic_Status",
    group = "R_T_N"
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

is_forbidden_comparison <- function(pair) {
  any(vapply(
    NOT_comparisons,
    function(x) identical(as.character(pair), as.character(x)),
    logical(1)
  ))
}

# -----------------------------
# 3. Plot one cell type
# -----------------------------
plot_one_type <- function(i) {
  
  subdata <- group_df_long %>%
    dplyr::filter(variable == total_types[i]) %>%
    dplyr::mutate(
      Group_Label = factor(Var3, levels = orders2)
    )
  
  maxx <- max(subdata$value, na.rm = TRUE)
  
  onetype_plot <- ggplot(
    subdata,
    aes(x = Group_Label, y = value, fill = Group_Label)
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
    geom_jitter(color = "black", size = 0.4, alpha = 0.9,height = 0) +
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
      
      if (!is_forbidden_comparison(pair)) {
        
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

s3q <- patchwork::wrap_plots(
  plot_list_with_legend,
  ncol = 3
) +
  patchwork::plot_layout(guides = "keep") +
  patchwork::plot_annotation() &
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5),
    panel.spacing = unit(0.1, "lines")
  )


ggsave(filename = "s3q.png", plot = s3q, width = 6, height = 4, units = "in")


