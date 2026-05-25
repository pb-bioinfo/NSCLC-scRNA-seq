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

setwd("~/diskD/lung/paneltotal/fig2/")

obj_epi <-qread("obj_epi_normal.qs")
obj_mali <-qread("obj_mali.qs")



color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')



#####

f2a <-  DimPlot(obj_epi, reduction = "umap",split.by = "R_T_N" , label=F,group.by ="Epi_SciAdv2020",ncol = 2)+labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = color.liberal)+theme(strip.text = element_text(size = 15))
ggsave(filename = "f2a.tiff", plot = f2a, width = 7, height = 3, units = "in")


#####

pt <- table(obj_epi@meta.data$R_T_N, obj_epi@meta.data$Epi_SciAdv2020)
pt <- as.data.frame(pt)




f2b1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))




f2b2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=14))






combined_plot <- (f2b1 | f2b2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot

ggsave(filename = "f2b.tiff", plot = combined_plot, width = 4.5, height = 4, units = "in")


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


meta_df <- obj_epi@meta.data %>%
  dplyr::select(
    Sample_ID,
    celltype = "Epi_SciAdv2020",
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


f2c <- final_plot

ggsave(filename = "f2c.tiff", plot = f2c, width = 5.1, height = 6, units = "in")



#####
Idents(obj_epi) = "Race_TNM"


f2d <- DEenrichRPlot(
  obj_epi,
  ident.1 = "B_N",
  ident.2 = "W_N",
  balanced = TRUE,
  logfc.threshold = 0.1, # Lowered slightly to catch more genes for pathway analysis
  assay = "RNA",         # Usually safer to specify "RNA"
  max.genes = 500,       # Increased: KEGG needs more genes than Hallmark to trigger hits
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "MSigDB_Hallmark_2020",
  num.pathway = 10,      # Increased to see more results
  return.gene.list = FALSE)


f2d[[1]] <- f2d[[1]]  + 
  labs(title =paste0("Black Normal ","Epithelial"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


f2d[[2]] <- f2d[[2]]  + 
  labs(title =paste0("White Normal ","Epithelial"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))



ggsave(filename = "f2d.tiff", plot = f2d, width = 6.2, height = 3, units = "in")

#####

Idents(obj_mali) = "Race_TNM"

f2f <- DEenrichRPlot(
  obj_mali,
  ident.1 = "B_T",
  ident.2 = "W_T",
  balanced = TRUE,
  logfc.threshold = 0.1, 
  assay = "RNA",         
  max.genes = 500,      
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "MSigDB_Hallmark_2020",
  num.pathway = 10,      
  return.gene.list = FALSE)
f2f


f2f[[1]] <- f2f[[1]]  + 
  labs(title =paste0("Black Tumor ","Malignant"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))



f2f[[2]] <- f2f[[2]]  + 
  labs(title =paste0("White Tumor ","Malignant"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))




ggsave(filename = "f2f.tiff", plot = f2f, width = 6.2, height = 3, units = "in")



