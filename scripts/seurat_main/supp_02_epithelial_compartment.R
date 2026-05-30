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
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


setwd("~/diskD/lung/paneltotal/fig2/")

obj_epi <-qread("obj_epi_normal.qs")
obj_mali <-qread("obj_mali.qs")
obj_epi_T <-qread("obj_epi_T.qs")



#####


s2a <- DimPlot(obj_epi, reduction = "umap", label=T, cols = col_vector,group.by = "RNA_snn_res.0.5",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")
s2a
ggsave(filename = "s2a.png", plot = s2a, width = 3.5, height = 3, units = "in")


#####



shortlist<- c("EPCAM", "AGER","SFTPA1","SCGB1A1","KRT17","TPPP3","SPI1","SPINK1")

plots <- lapply(1:8, function(i) {
  if (i == 5){
    FeaturePlot(obj_epi,features = shortlist[i], reduction = "umap",cols = c('grey','red'))+labs(y = "UMAP_2", x = "UMAP_1",max.cutoff = 4.6)+theme(legend.position = "none")
  }else{
    FeaturePlot(obj_epi,features = shortlist[i], reduction = "umap",cols = c('grey','red'))+labs(y = " ", x = " ",max.cutoff = 4.6,no.axes = TRUE)+theme(legend.position = "none")
  }
})
s2b <- wrap_plots(plots, ncol = 4) +  plot_layout(guides = 'collect') & 
  guides(fill = guide_legend(title = "Expression"))
#  theme(legend.position = "right") &

ggsave(filename = "s2b.png", plot = s2b, width = 7.2, height = 3.6,units = "in")


#####

s2c <- DimPlot(obj_epi, reduction = "umap", label=T, cols = color.liberal,group.by = "Epi_SciAdv2020",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")
s2c

ggsave(filename = "s2c.png", plot = s2c, width = 3.5, height = 3, units = "in")


#####


#####


pt <- table(obj_epi@meta.data$Sample_ID, obj_epi@meta.data$Epi_SciAdv2020)
pt <- as.data.frame(pt)




s2d1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))




s2d2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))






combined_plot <- (s2d1 | s2d2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot

ggsave(filename = "s2d.png", plot = combined_plot, width = 4.5, height = 4, units = "in")


#####


custom_colors <- c("orange",      # orange
                   "#377EB8",      # blue
                   "#A50F15",      # dark red
                   "#9ECAE1")


# -----------------------------
# 1. Build patient-level proportions
# -----------------------------
meta_df <- obj_epi@meta.data %>%
  dplyr::select(
    Sample_ID,
    celltype = "Epi_SciAdv2020",
    group = "Type2"
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



s2e <- final_plot

ggsave(filename = "s2e.png", plot = s2e, width = 6, height = 6, units = "in")


#####



pt <- table(obj_epi_T@meta.data$Sample_ID, obj_epi_T@meta.data$maligant_sciadv2020)
pt <- as.data.frame(pt)




s2g1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s2g2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))






combined_plot <- (s2g1 / s2g2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot


ggsave(filename = "s2g.png", plot = combined_plot, width = 9, height = 5, units = "in")

#####

s2h <- DimPlot(obj_mali, reduction = "umap", label=T, cols = col_vector,group.by = "RNA_snn_res.0.5",label.size = 4) +
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")+ NoLegend()
s2h


ggsave(filename = "s2h.png", plot = s2h, width = 4, height = 3.5, units = "in")


#####


s2i <- DimPlot(obj_mali, reduction = "umap", label=F, cols = col_vector,group.by = "Sample_ID") +
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")
s2i
ggsave(filename = "s2i.png", plot = s2i, width = 8, height = 4.5, units = "in")


#####

s2j <- DimPlot(obj_mali, reduction = "umap",split.by = "Type2" , label=F,group.by ="RNA_snn_res.0.5",ncol = 2)+ NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = col_vector)

ggsave(filename = "s2j.png", plot = s2j, width = 5.6, height = 3, units = "in")



#####


s2k <- DimPlot(obj_mali, reduction = "umap",split.by = "R_T_N" , label=F,group.by ="RNA_snn_res.0.5",ncol = 2)+ NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = col_vector)


ggsave(filename = "s2k.png", plot = s2k, width = 5.6, height = 3, units = "in")


