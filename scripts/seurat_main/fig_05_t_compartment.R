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

setwd("~/diskD/lung/paneltotal/fig5/")

obj_T <-qread("obj_T.qs")
obj_CD8 <-qread("obj_CD8.qs")
obj_NK <-qread("obj_NK.qs")


color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

#####

f5a <- DimPlot(obj_T, reduction = "umap.t.rpca",split.by = "R_T_N" , label=F,group.by ="stage2",ncol = 4)+labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = color.liberal)+theme(strip.text = element_text(size = 15))+ NoLegend()
f5a
ggsave(filename = "f5a.png", plot = f5a, width = 10, height = 3, units = "in")




#####



pt <- table(obj_T@meta.data$R_T_N, obj_T@meta.data$stage2v2)
pt <- as.data.frame(pt)



f5b1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,,size = 12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




f5b2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))






combined_plot <- (f5b1 | f5b2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
combined_plot
ggsave(filename = "f5b.png", plot = combined_plot, width = 6, height = 4.5, units = "in")


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
meta_df <- obj_T@meta.data %>%
  dplyr::select(
    Sample_ID,
    celltype = "stage2v2",
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



f5c <- final_plot

ggsave(filename = "f5c.png", plot = f5c, width = 4.5, height = 9, units = "in")

#####

Idents(obj_CD8) = "Race_TNM"

f5d <- DEenrichRPlot(
  obj_CD8,
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

f5d[[1]] <- f5d[[1]]  + 
  labs(title =paste0("Black Tumor ","CD8 T"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)   
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


f5d[[2]] <- f5d[[2]]  + 
  labs(title =paste0("White Tumor ","CD8 T"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)   
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


ggsave(filename = "f5d.png", plot = f5d, width = 6, height = 3, units = "in")



#####

#CD8 Exhausted 
obj_exhausted <-qread("obj_exhausted.qs")


Idents(obj_exhausted) = "Race_TNM"
f5e <- DEenrichRPlot(
  obj_exhausted,
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




f5e[[1]] <- f5e[[1]]  + 
  labs(title =paste0("Black Tumor ","CD8 Exhausted"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)   
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))



f5e[[2]] <- f5e[[2]] + 
  labs(title =paste0("White Tumor ","CD8 Exhausted"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)   
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


ggsave(filename = "f5e.png", plot = f5e, width = 6, height = 3, units = "in")


#####

#Signature Scores
cytotoxic2 = c("GNLY","IFNG","PRF1","GZMK","GZMH","GZMM","KLRK1", "KLRB1", "KLRD1", "CTSW","CST7")
exhausted_genes <- c("PDCD1","LAG3","TIGIT","HAVCR2","BTLA")

sce_score <- AddModuleScore(
  object = obj_CD8,
  features = list(cytotoxic2),
  name = 'cytotoxic_score'
)


sce_score <- AddModuleScore(
  object = sce_score,
  features = list(exhausted_genes),
  name = 'exhausted_score'
)



metadata2<- sce_score@meta.data[,c("R_T_N","stage2", "cytotoxic_score1","exhausted_score1","Race_TNM")]

f5f1<- ggplot(metadata2, aes(x = cytotoxic_score1, colour = R_T_N)) +scale_color_manual(values = color.liberal)+
  stat_ecdf(geom = "step", size = 0.5) + 
  labs(
    title = "",
    x = "Cytotoxic Score",
    y = "Cumulative Proportion",
    colour = " "
  )+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size =10),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 10))+
  theme(axis.line = element_line(colour = "black"))

f5f2<- ggplot(metadata2, aes(x = exhausted_score1, colour = R_T_N)) +scale_color_manual(values = color.liberal)+
  stat_ecdf(geom = "step", size = 0.5) +
  labs(
    title = "",
    x = "Exhaustion Score",
    y = "",
    colour = " "
  )+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5,size =10), 
        panel.background = element_blank(),
        axis.text.y = element_text(size = 10))+
  theme(axis.line = element_line(colour = "black"))



combined_plot <- (f5f1 | f5f2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")


ggsave(filename = "f5f.png", plot = combined_plot, width = 7.5, height = 3, units = "in")


ks.test(metadata2$cytotoxic_score1[metadata2$Race_TNM == "B_T"], metadata2$cytotoxic_score1[metadata2$Race_TNM == "W_T"])

ks.test(metadata2$exhausted_score1[metadata2$Race_TNM == "B_T"], metadata2$exhausted_score1[metadata2$Race_TNM == "W_T"])

#####



state_time_ver2_total <- qread("~/diskD/lung/Tmonocole.qs")

cell_tensor = table(state_time_ver2_total$stage2,state_time_ver2_total$State,state_time_ver2_total$Race_TNM)


exhausted_T_df <- as.data.frame.matrix(cell_tensor["CD8 Exhausted", , c("B_T", "W_T")]) %>%
  rownames_to_column(var = "State") %>%
  pivot_longer(
    cols = c("B_T", "W_T"),
    names_to = "Race_TNM",
    values_to = "Cell_Count"
  )


exhausted_proportions_df <- exhausted_T_df %>%
  dplyr::group_by(Race_TNM) %>% # Group by B_T and W_T
  dplyr::mutate(Proportion = Cell_Count / sum(Cell_Count)) %>% # Calculate proportion within each group
  dplyr::ungroup()  



cytotoxic_T_df <- as.data.frame.matrix(cell_tensor["CD8 Cytotoxic", , c("B_T", "W_T")]) %>%
  rownames_to_column(var = "State") %>%
  pivot_longer(
    cols = c("B_T", "W_T"),
    names_to = "Race_TNM",
    values_to = "Cell_Count"
  )


cytotoxic_proportions_df <-cytotoxic_T_df %>%
  dplyr::group_by(Race_TNM) %>% # Group by B_T and W_T
  dplyr::mutate(Proportion = Cell_Count / sum(Cell_Count)) %>% # Calculate proportion within each group
  dplyr::ungroup()  



activated_T_df <- as.data.frame.matrix(cell_tensor["CD8 Activated", , c("B_T", "W_T")]) %>%
  rownames_to_column(var = "State") %>%
  pivot_longer(
    cols = c("B_T", "W_T"),
    names_to = "Race_TNM",
    values_to = "Cell_Count"
  )

activated_proportions_df <- activated_T_df %>%
  dplyr::group_by(Race_TNM) %>% # Group by B_T and W_T
  dplyr::mutate(Proportion = Cell_Count / sum(Cell_Count)) %>% # Calculate proportion within each group
  dplyr::ungroup()  


proliferating_T_df <- as.data.frame.matrix(cell_tensor["CD8 Proliferating", , c("B_T", "W_T")]) %>%
  rownames_to_column(var = "State") %>%
  pivot_longer(
    cols = c("B_T", "W_T"),
    names_to = "Race_TNM",
    values_to = "Cell_Count"
  )

proliferating_proportions_df <- proliferating_T_df %>%
  dplyr::group_by(Race_TNM) %>% # Group by B_T and W_T
  dplyr::mutate(Proportion = Cell_Count / sum(Cell_Count)) %>% # Calculate proportion within each group
  dplyr::ungroup()  


p_exhausted2 <- ggplot(exhausted_proportions_df, aes(x = State, y = Proportion, fill = Race_TNM)) +
  geom_col(position = "dodge") +  
  scale_fill_manual(values = c("B_T" = "#1f77b4", "W_T" = "#ff7f0e")) +   
  labs(
    title = "CD8+ T Exhausted",    
    x = NULL,     
    y = "Proportion of Cells",
    fill = "Patient Group"
  ) +
  theme_classic() +  
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))+
  theme(axis.text.x = element_text(size = 12))+
  theme(
    legend.position = c(0.7, 0.95), # Position legend inside the plot (x=0.95, y=0.95)
    legend.justification = c("right", "top"), # Anchor the legend's top-right corner
    legend.background = element_rect(color = "black", linewidth = 0.5) # Draw a box
  )




p_cytotoxic2 <- ggplot(cytotoxic_proportions_df, aes(x = State, y = Proportion, fill = Race_TNM)) +
  geom_col(position = "dodge") +  
  scale_fill_manual(values = c("B_T" = "#1f77b4", "W_T" = "#ff7f0e")) +   
  labs(
    title = "CD8+ T Cytotoxic",    
    x = NULL,     
    y = "Proportion of Cells",
    fill = "Patient Group"
  ) +
  theme_classic() +  
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(size = 12))+
  theme( panel.background = element_blank(),text = element_text(size=12))+
  NoLegend()


p_activated2 <- ggplot(activated_proportions_df, aes(x = State, y = Proportion, fill = Race_TNM)) +
  geom_col(position = "dodge") +  
  scale_fill_manual(values = c("B_T" = "#1f77b4", "W_T" = "#ff7f0e")) +   
  labs(title = "CD8+ T Activated", x = NULL, y = "Proportion of Cells") +
  theme_classic() +  
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(size = 12))+
  theme( panel.background = element_blank(),text = element_text(size=12))+
  NoLegend()



p_proliferating2 <- ggplot(proliferating_proportions_df, aes(x = State, y = Proportion, fill = Race_TNM)) +
  geom_col(position = "dodge") +  
  scale_fill_manual(values = c("B_T" = "#1f77b4", "W_T" = "#ff7f0e")) +   
  labs(title = "CD8+ T Proliferating", x = "Monocle State", y = "Proportion of Cells") +
  theme_classic() +  
  theme_bw(base_size = 12) +
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(size = 12))+
  theme( panel.background = element_blank(),text = element_text(size=12))+
  NoLegend()
# --- Combine the three plots vertically and create a single legend ---
final_plot2 <- p_exhausted2 + p_cytotoxic2 + p_activated2 + p_proliferating2+
  plot_layout(ncol = 2) 

f5h <- final_plot2 


ggsave(filename = "f5h.png", plot = f5h,width = 10, height = 4, units = "in")


#####

f5i <- DimPlot(obj_NK, reduction = "umap.NK.rpca", label=T,label.size = 4, repel= TRUE,cols = color.liberal,group.by = "NK_subtype") + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")
f5i


ggsave(filename = "f5i.png", plot = f5i, width = 3.8, height = 3.5, units = "in")


#####




sce_score <- AddModuleScore(
  object = obj_NK,
  features = list(cytotoxic2),
  name = 'cytotoxic_score'
)



metadata2<- sce_score@meta.data[,c("Race_Type","NK_subtype","cytotoxic_score1")]

f5j<- ggplot(metadata2, aes(x = cytotoxic_score1, colour = Race_Type)) +scale_color_manual(values = color.liberal)+
  stat_ecdf(geom = "step", size = 0.5) + # Use stat_ecdf with geom = "step"
  labs(
    title = "",
    x = "Cytotoxic Score",
    y = "Cumulative Proportion",
    colour = " "
  )+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, size =10),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 10))+
  theme(axis.line = element_line(colour = "black"))


ggsave(filename = "f5j.png", plot = f5j, width = 5.6, height = 3, units = "in")



ks.test(metadata2$cytotoxic_score1[metadata2$Race_Type == "B_NK_CD16_positive"], metadata2$cytotoxic_score1[metadata2$Race_Type == "W_NK_CD16_positive"])

ks.test(metadata2$cytotoxic_score1[metadata2$Race_Type == "B_NK_CD16_intermediate"], metadata2$cytotoxic_score1[metadata2$Race_Type == "W_NK_CD16_intermediate"])

ks.test(metadata2$cytotoxic_score1[metadata2$Race_Type == "B_NK_CD16_negative"], metadata2$cytotoxic_score1[metadata2$Race_Type == "W_NK_CD16_negative"])

