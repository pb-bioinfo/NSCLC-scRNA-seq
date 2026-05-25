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

#####

f4a <- DimPlot(obj_myeloid, reduction = "umap.Myeloid.rpca",split.by = "R_T_N" , label=F,group.by ="stage1_v3",ncol = 2)+labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = color.liberal)+theme(strip.text = element_text(size = 15))+NoLegend()

ggsave(filename = "f4a.tiff", plot = f4a, width = 6, height = 6, units = "in")

#####





pt <- table(obj_myeloid@meta.data$R_T_N, obj_myeloid@meta.data$stage1_v3)
pt <- as.data.frame(pt)



f4b1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




f4b2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))


combined_plot <- (f4b1 | f4b2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
combined_plot

ggsave(filename = "f4b.tiff", plot = combined_plot, width = 5, height = 4, units = "in")



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
meta_df <- obj_myeloid@meta.data %>%
  dplyr::select(
    Sample_ID,
    celltype = "stage1_v3",
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


f4c <- final_plot

ggsave(filename = "f4c.tiff", plot = f4c, width = 5.1, height = 4.5, units = "in")


#####

Idents(obj_macro)<-"Race_TNM"
f4d <- DEenrichRPlot(
  obj_macro,
  ident.1 = "B_T",
  ident.2 = "W_T",
  balanced = TRUE,
  logfc.threshold = 0.1,
  assay = "RNA",         
  max.genes = 500,       
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "KEGG_2021_Human", 
  num.pathway = 10,      # Increased to see more results
  return.gene.list = FALSE
)


f4d[[1]] <- f4d[[1]]  + 
  labs(title =paste0("Black Tumor ","TAM"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


f4d[[2]] <- f4d[[2]]  + 
  labs(title =paste0("White Tumor ","TAM"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))



ggsave(filename = "f4d.tiff", plot = f4d, width = 7.1, height = 3, units = "in")



#####


M1_gene = c("IFNG","TNF", "NOS2","IDO1","CALHM6","CXCL9","CD86","IL12A","IL23A")

M2_gene = c("IL4","IL13", "MRC1","CD209","TGM2","CCL24","CCL26","ARG1", "TGFBR1", "TGFB1", "PPARG", "CHI3L1", "STAB1", "MSR1","ALDH1A2","SPP1")

sce_score <- AddModuleScore(
  object = obj_macro,
  features = list(M1_gene),
  name = 'M1_score'
)
sce_score <- AddModuleScore(
  object = sce_score,
  features = list(M2_gene),
  name = 'M2_score'
)



metadata2<- sce_score@meta.data[,c("R_T_N","M1_score1","M2_score1")]

f4e1<- ggplot(metadata2, aes(x = M1_score1, colour = R_T_N)) +scale_color_manual(values = color.liberal)+
  stat_ecdf(geom = "step", size = 0.5) + # Use stat_ecdf with geom = "step"
  labs(
    title = "",
    x = "M1 Score",
    y = "Cumulative Proportion",
    colour = " "
  )+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5,size =10), 
        panel.background = element_blank(),
        axis.text.y = element_text(size = 10))+
  theme(axis.line = element_line(colour = "black"))



f4e2<- ggplot(metadata2, aes(x = M2_score1, colour = R_T_N)) +scale_color_manual(values = color.liberal)+
  stat_ecdf(geom = "step", size = 0.5) + # Use stat_ecdf with geom = "step"
  labs(
    title = "",
    x = "M2 Score",
    y = "",
    colour = " "
  )+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5,size =10), 
        panel.background = element_blank(),
        axis.text.y = element_text(size = 10))+
  theme(axis.line = element_line(colour = "black"))




combined_plot <- (f4e1 | f4e2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot

ggsave(filename = "f4e.tiff", plot = combined_plot,width = 6.3, height = 2.5, units = "in")

ks.test(metadata2$M1_score1[metadata2$Race_TNM == "B_T"], metadata2$M1_score1[metadata2$Race_TNM == "W_T"])
ks.test(metadata2$M2_score1[metadata2$Race_TNM == "B_T"], metadata2$M2_score1[metadata2$Race_TNM == "W_T"])


#####
antigen_presentation_genes <- c(
"CD1B", "CD1D", "CD40", "CD80", "CD86",
"HLA-A", "HLA-B", "HLA-C",
"TAP1", "TAP2", "ERAP1", "ERAP2",   
"HLA-DRA", "HLA-DRB1", "HLA-DMA", "HLA-DMB", "HLA-DQA1", "HLA-DQB1"
)


Idents(obj_macro) <-  "Race_TNM"

# Here we run standard DE to get adjusted p-values properly.
markers <- FindMarkers(
  obj_macro,
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
f4f <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw() +
  labs(
    title = "Enriched in black tumor TAMs",
    x = "",
    y = "Average Log2FC",
    size = "-log10(p-adj)",
    color = "Log2FC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.box = "horizontal")+
  guides( color = "none")

ggsave(filename = "f4f.tiff", plot = f4f, width = 4.5, height = 3.5, units = "in")



#####

myeloid_tensor = table(trajectory_subset$stage1_v3,trajectory_subset$State,trajectory_subset$Race_TNM)
unique_types =  c("Monocyte","Mono-mac","TAM")

plot_list2 <- list()

for (type in unique_types) {

  cell_type_name <- type
  
  cell_type_slice <- myeloid_tensor[cell_type_name, , c("B_T", "W_T")]
  
  cell_type_df <- as.data.frame.matrix(cell_type_slice) %>%
    rownames_to_column(var = "State") %>%
    pivot_longer(
      cols = c("B_T", "W_T"),
      names_to = "Race_TNM",
      values_to = "Cell_Count"
    )
  cell_proportion_df <- cell_type_df %>%
    group_by(Race_TNM) %>% # Group by B_T and W_T
    mutate(Proportion = Cell_Count / sum(Cell_Count)) %>% # Calculate proportion within each group
    ungroup() 
  
  p <- ggplot(cell_proportion_df, aes(x = State, y = Proportion, fill = Race_TNM)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("B_T" = "#1f77b4", "W_T" = "#ff7f0e")) +
    labs(
      title = cell_type_name, # Use the cell type name as the title
      x = NULL,
      y = "Proportion of Cells",
      fill = "Patient Group"
    ) +
    theme_classic() +
    theme_bw(base_size = 15) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.background = element_blank(), text = element_text(size = 12))+
    theme(axis.text.x = element_text(size = 12))

  plot_list2[[type]] <- p
  
}

final_stacked_plot2 <- wrap_plots(plot_list2, ncol = 1) +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'right') 


ggsave(filename = "f4h.tiff", plot = final_stacked_plot2, width = 6, height = 6.5, units = "in")



#####


f4i <- DimPlot(obj_DC, reduction = "umap", label=T, cols = color.liberal,group.by = "DC",split.by = "R_T_N" ,
               ncol = 4,label.size = 4,repel = TRUE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+labs(title = "")+
  theme(
    strip.text = element_text(size = 15, face = "bold")
  )

ggsave(filename = "f4i.tiff", plot = tiff, width = 10, height = 3.5, units = "in")
