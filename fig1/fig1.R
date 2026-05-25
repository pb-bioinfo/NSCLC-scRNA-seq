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
library(purrr)

library(data.table)
library(reshape2)
library(viridis)
library(ggsci)
library(qusage)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(DescTools)

library(pheatmap)
library(stats)
library(lsa)
library(grid)
library(dplyr)


library(Seurat)
library(igraph)

#library(plyr)

library(msigdbr)
library(presto)
library(MAST)
library(openxlsx)
library(CellChat)
setwd("~/diskD/lung/paneltotal/fig1/")


rds1 <-qread("rds1.qs")

color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')


##### Figure 1A #####


f1a <- DimPlot(rds1, reduction = "umap.rpca",split.by = "R_T_N" , label=F,group.by ="celltype2_RNA_snn_res.2",ncol = 4)+labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = color.liberal)+theme(strip.text = element_text(size = 15))

ggsave(filename = "f1a.tiff", plot = f1a, width = 11, height = 3, units = "in")



#####


pt <- table(rds1@meta.data$R_T_N, rds1@meta.data$celltype2_RNA_snn_res.2)
pt <- as.data.frame(pt)




f1b1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))




f1b2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))






combined_plot <- (f1b1 | f1b2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot
ggsave(filename = "f1b.tiff", plot = combined_plot, width = 5.5, height = 4, units = "in")



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


meta_df <- rds1@meta.data %>%
  dplyr::select(
    Sample_ID,
    celltype = "celltype2_RNA_snn_res.2",
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


final_plot



f1c <- final_plot

ggsave(filename = "f1c.tiff", plot = f1c, width = 4.5, height = 6, units = "in")


#rm(f1a,f1b,f1c)