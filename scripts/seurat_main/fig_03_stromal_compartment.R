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

library(dplyr)
library(Seurat)
library(presto)
library(enrichR)

color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')


setwd("~/diskD/lung/paneltotal/fig3/")

obj_Fibro <-qread("obj_Fibro.qs")

obj_caf2 <- qread("obj_caf2.qs")

obj_Endo <- qread("obj_Endo.qs")
#####


f3a <- DimPlot(obj_Fibro, reduction = "umap.endo.rpca",split.by = "R_T_N" , label=F,group.by ="endo",ncol = 4)+labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = color.liberal)
ggsave(filename = "f3a.png", plot = f3a, width = 11, height = 3, units = "in")




#####




pt <- table(obj_Fibro@meta.data$R_T_N, obj_Fibro@meta.data$endo)
pt <- as.data.frame(pt)



f3b1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))




f3b2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))





combined_plot <- (f3b1 | f3b2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
combined_plot
ggsave(filename = "f3b.png", plot = combined_plot, width = 4.8, height = 4, units = "in")



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
meta_df <- obj_Fibro@meta.data %>%
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

f3c <- patchwork::wrap_plots(
  plot_list_with_legend,
  ncol = 3
) +
  patchwork::plot_layout(guides = "keep") +
  patchwork::plot_annotation() &
  theme(
    plot.margin = margin(0.5, 0.5, 0.5, 0.5),
    panel.spacing = unit(0.1, "lines")
  )


ggsave(filename = "f3c.png", plot = f3c, width = 6, height = 5, units = "in")


#####


Idents(obj_caf2) <- "Race_TNM"


f3d <- DEenrichRPlot(
  obj_caf2,
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
  return.gene.list = FALSE)



f3d[[1]] <- f3d[[1]] + 
  labs(title =paste0("Black Tumor ","CAF-2"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.9) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )


f3d[[2]] <- f3d[[2]] + 
  labs(title =paste0("White Tumor ","CAF-2"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.9) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )

ggsave(filename = "f3d.png", plot = f3d, width = 7.5, height = 3, units = "in")




#####


# Lipid antigen presentation
lipid_antigen_genes <- c(
  "CD1D", "CD1B", "B2M", "PSAP"
)

# Co-stimulatory molecules
costimulatory_genes <- c(
  "CD40", "CD80", "CD86"
)

# MHC class I antigen presentation
mhc_class1_genes <- c(
  "HLA-A", "HLA-B", "HLA-C",
  "TAP1", "TAP2",
  "ERAP1", "ERAP2"
)

# MHC class II antigen presentation
mhc_class2_genes <- c(
  "HLA-DRA", "HLA-DRB1",
  "HLA-DMA", "HLA-DMB",
  "HLA-DQA1", "HLA-DQB1"
)

all_antigen_presentation_genes <- c(
  lipid_antigen_genes,
  costimulatory_genes,
  mhc_class1_genes,
  mhc_class2_genes
)

Idents(obj_caf2) <- "Race_TNM"

# Here we run standard DE to get adjusted p-values properly.
markers <- FindMarkers(
  obj_caf2,
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
f3e <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
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

ggsave(filename = "f3e.png", plot = f3e, width = 5, height = 3.5, units = "in")


#####


f3f <- DimPlot(obj_Endo, reduction = "umap.endo.rpca",split.by = "R_T_N" , label=F,group.by ="endo",ncol = 4)+labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = color.liberal)


ggsave(filename = "f3f.png", plot = f3f, width = 11, height = 3, units = "in")



#####



pt <- table(obj_Endo@meta.data$R_T_N, obj_Endo@meta.data$endo)
pt <- as.data.frame(pt)



f3g1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))




f3g2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))





combined_plot <- (f3g1 | f3g2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
combined_plot
ggsave(filename = "f3g.png", plot = combined_plot, width = 4.8, height = 4, units = "in")





#####


obj_TEC <- subset(obj_Endo,endo == "Tumor EC")

Idents(obj_TEC) <- "Race_TNM"

f3h <- DEenrichRPlot(
  obj_TEC,
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
  return.gene.list = FALSE)


f3h[[1]] <- f3h[[1]] + 
  labs(title =paste0("Black Tumor ","Tumor EC"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.9) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )


f3h[[2]] <- f3h[[2]] + 
  labs(title =paste0("White Tumor ","Tumor EC"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.9) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )
ggsave(filename = "f3h.png", plot = f3h, width = 7.3, height = 3, units = "in")



