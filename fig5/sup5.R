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

#library(plyr)

library(presto)
library(MAST)
library(CellChat)

setwd("~/diskD/lung/paneltotal/fig3/")

obj_T <-qread("obj_T.qs")
obj_CD8 <-qread("obj_CD8.qs")

obj_NK <-qread("obj_NK.qs")



color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


s5a <- DimPlot(obj_T, reduction = "umap.t.rpca", label=T, cols = col_vector,group.by = "seurat_clusters",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")
s5a
ggsave(filename = "s5a.tiff", plot = s5a, width = 4.5, height = 4, units = "in")




shortlist<- c("CD3D", "CD8A","CD4","FOXP3","PRF1","TIGIT","MKI67","TRGC1","CCR7")


plots <- lapply(1:9, function(i) {
  if (i == 6){
    FeaturePlot(obj_T,features = shortlist[i], reduction = "umap.t.rpca",cols = c('grey','red'))+labs(y = "UMAP_2", x = "UMAP_1",max.cutoff = 4.6)+theme(legend.position = "none")
  }else{
    FeaturePlot(obj_T,features = shortlist[i], reduction = "umap.t.rpca",cols = c('grey','red'))+labs(y = " ", x = " ",max.cutoff = 4.6,no.axes = TRUE)+theme(legend.position = "none")
  }
})
s5b <- wrap_plots(plots, ncol = 5) +  plot_layout(guides = 'collect') +
  guides(fill = guide_legend(title = "Expression")) +
  theme(legend.position = "right")  

ggsave(filename = "s5b.tiff", plot = s5b, width = 8.2, height = 3.6,units = "in")




#####

avgexpT2 <- AverageExpression(obj_T, return.seurat = T, group.by = 'seurat_clusters')

myorders <- c("g6", "g10","g12",
"g21", "g11", "g44",
"g3","g7", "g33", "g38",
"g46", "g20",
"g2","g5",
"g29", "g43", "g42", "g9", "g8", "g13" ,
"g19",
"g14","g37",
"g24",  "g15",  "g0", "g17", "g36", "g1", "g25", "g28", "g40", "g18",
"g31", "g23", "g39", "g22", "g34",
"g30", "g27",
"g26", "g16", "g32", "g41", "g45",
"g4", "g35"
)

general_T = c("CD3G",'CD3D',"CD3E", "CD4", "CD8A","CD8B")
Gamma_delta = c( "TRDC","TRGC1","TRGC2")
act1 = c("CD40LG","CD69","CD38","HLA-DRA","HLA-DRB1","ICOS","TNFRSF4","TNFRSF9")
exhau_inhib = c("PDCD1","LAG3","TIGIT","HAVCR2","BTLA")
cytotoxic2 = c("GNLY","IFNG","PRF1","GZMK","GZMH","GZMM","KLRK1", "KLRB1", "KLRD1", "CTSW","CST7")
tfh = c("CXCL13","MAF","CXCR5")
th1 = c("CXCR3","STAT4","RGS1","IL12RB2")
th2 = c('GATA3',"STAT6","IL4","CCR3","CCR4","CCR8")
th9 = c("IL9")
th17 = c("IRF4","CCR6","IL17A","IL17F")
th22 = c("IL22")
Treg = c('FOXP3',  'CTLA4',"IL2RA", 'TNFRSF18', 'IKZF2')
Prolif = c("MKI67","CDK1","STMN1")
naive_memory2 = c( "SELL","IL7R", "CCR7","LEF1","TCF7","IL2RA")
markers2 = c(general_T,act1,tfh,th1,th17,th22,Treg,th9,cytotoxic2,exhau_inhib,Prolif,Gamma_delta,naive_memory2)



Idents(avgexpT2) <- factor(Idents(avgexpT2), levels = myorders)

s5c <- DoHeatmap(
  object = avgexpT2,
  features = markers2,  # Genes
  size = 5,
  draw.lines = FALSE,
  disp.min = 0,
  disp.max = 2,
  group.by = "seurat_clusters"  
) +
  scale_fill_gradientn(colors = c("lightgrey", "blue")) +
  scale_x_discrete(labels = gsub("^g", "", levels(avgexpT2$seurat_clusters_clean))) +  # Explicitly remove "g" prefix
  theme(
    axis.text.y = element_text(size = 10, family = "Arial", color = "black"),  # Gene labels
    axis.text.x.bottom = element_text(size = 8, angle = 90,hjust = 1, family = "Arial", color = "black"),  # Bottom axis cluster labels
    axis.text.x.top = element_blank(),  #axis.text.x.top = element_text(size = 8, angle = 45,  family = "Arial",hjust = 1, family = "Arial", color = "black"),  # Top axis cluster labels
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank()   # Remove y-axis title
  )



ggsave(filename = "s5c.tiff", plot = s5c, width = 6.5, height = 9, units = "in")


#####


s5d <- DimPlot(obj_T, reduction = "umap.t.rpca", label=F, cols = color.liberal,group.by = "stage2v2",repel = FALSE)+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")+
  theme(legend.text = element_text(size = 12))
s5d
ggsave(filename = "s5d.tiff", plot = s5d, width = 4.5, height = 3, units = "in")


#####


pt <- table(obj_T@meta.data$Sample_ID, obj_T@meta.data$stage2v2)
pt <- as.data.frame(pt)




s5e1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x =element_blank())+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s5e2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))



combined_plot <- (s5e1 / s5e2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")


ggsave(filename = "s5e.tiff", plot = combined_plot, width = 11, height = 5, units = "in")



#####
s5f <- DimPlot(obj_T, reduction = "umap.t.rpca",split.by = "T_N" , label=F,group.by ="stage2v2",ncol = 2)+ NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = color.liberal)+
  theme(strip.text = element_text(face = "bold", size = 15))

ggsave(filename = "s5f.tiff", plot =s5f, width = 5.6, height = 3, units = "in")


#####




pt <- table(obj_T@meta.data$T_N, obj_T@meta.data$stage2v2)
pt <- as.data.frame(pt)



s5g1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,,size = 12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s5g2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 12) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 12))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))






combined_plot <- (s5g1 | s5g2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
combined_plot
ggsave(filename = "s5g.tiff", plot = combined_plot, width = 6, height = 4.5, units = "in")


#####




s5h <- DimPlot(obj_T, reduction = "umap.t.rpca",split.by = "Race2" , label=F,group.by ="stage2v2",ncol = 2)+ NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+
  ggtitle(NULL)+scale_color_manual(values = color.liberal)+
  theme(strip.text = element_text(face = "bold", size = 15))

ggsave(filename = "s5h.tiff", plot = s5h, width = 5.6, height = 3, units = "in")


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
    Sample_ID,
    celltype = "stage2v2",
    group = "Race2"
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


s5i <- final_plot
ggsave(filename = "s5i.tiff", plot = s5i, width = 7, height = 8, units = "in")

#####
Idents(obj_CD8) <- "Race_TNM"
s5j <- DEenrichRPlot(
  obj_CD8,
  ident.1 = "B_T",
  ident.2 = "W_T",
  balanced = TRUE,
  logfc.threshold = 0.1,  
  assay = "RNA",           
  max.genes = 500,         
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "KEGG_2021_Human", # <--- THIS IS THE CHANGE
  num.pathway = 10,        
  return.gene.list = FALSE
)

s5j[[1]] <- s5j[[1]]  + 
  labs(title =paste0("Black Tumor ","CD8 T"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


s5j[[2]] <- s5j[[2]]  + 
  labs(title =paste0("White Tumor ","CD8 T"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))

ggsave(filename = "s5j.tiff", plot = s5j, width =  7.5, height = 3, units = "in")


#####


Idents(obj_exhausted)<- "Race_TNM"
s5k <- DEenrichRPlot(
  obj_exhausted,
  ident.1 = "B_T",
  ident.2 = "W_T",
  balanced = TRUE,
  logfc.threshold = 0.1,  
  assay = "RNA",           
  max.genes = 500,         
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "KEGG_2021_Human", # <--- THIS IS THE CHANGE
  num.pathway = 10,        
  return.gene.list = FALSE
)

s5k[[1]] <- s5k[[1]]  + 
  labs(title =paste0("Black Tumor ","CD8 Exhausted"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))

s5k[[2]] <- s5k[[2]]  + 
  labs(title =paste0("White Tumor ","CD8 Exhausted"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


ggsave(filename = "s5k.tiff", plot = s5k, width = 7.5, height = 3, units = "in")


#####


pt <- table(obj_CD8@meta.data$State, obj_CD8@meta.data$stage2)
pt <- as.data.frame(pt)




s5m1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s5m2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))







combined_plot <- (s5m1 | s5m2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot


ggsave(filename = "s5m.tiff", plot = combined_plot, width = 5.5, height = 3, units = "in")



#####



pt <- table(obj_CD8@meta.data$State, obj_CD8@meta.data$Race_TNM)
pt <- as.data.frame(pt)




s5n1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




s5n2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))







combined_plot <- (s5n1 | s5n2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot
ggsave(filename = "s5n.tiff", plot = combined_plot, width = 4.5, height = 3, units = "in")


#####





Idents(obj_exhausted)<- "State"

s5o <- DEenrichRPlot(
  obj_exhausted,
  ident.1 = "2",
  ident.2 = "5",
  balanced = TRUE,
  logfc.threshold = 0.25,   
  assay = "RNA",           
  max.genes = 100,         
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "MSigDB_Hallmark_2020",
  num.pathway = 10,        
  return.gene.list = FALSE)




s5o[[1]] <- s5o[[1]]  + 
  labs(title =paste0("State 2 ","CD8 Exhausted"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14) 
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


s5o[[2]] <- s5o[[2]]  + 
  labs(title =paste0("State 5 ","CD8 Exhausted"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))



ggsave(filename = "s5o.tiff", plot = s5o, width = 6, height = 3, units = "in")



#####
obj_Treg <-qread("obj_Treg.qs")


Idents(obj_Treg)<- "State"

s5q <- DEenrichRPlot(
  obj_Treg,
  ident.1 = "4",
  ident.2 = "5",
  balanced = TRUE,
  logfc.threshold = 0.1,   
  assay = "RNA",           
  max.genes = 500,         
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  enrich.database = "KEGG_2021_Human", 
  num.pathway = 10,        
  return.gene.list = FALSE
)

s5q[[1]] <- s5q[[1]]  + 
  labs(title =paste0("State 4 ","Treg"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))


s5q[[2]] <- s5q[[2]]  + 
  labs(title =paste0("State 5 ","Treg"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.85) +theme(
    plot.title  = element_text(size = 14)  
  )+  theme(axis.text.x = element_text(size = 12),axis.title.x  =element_text(size=12))



ggsave(filename = "s5q.tiff", plot = s5q, width = 6.5, height = 3, units = "in")




#####



s5s = FeaturePlot(obj_NK,features = "FCGR3A", reduction = "umap.NK.rpca",cols = c('grey','red'))+
  labs(y = "UMAP_2", x = "UMAP_1",max.cutoff = 4.6)

ggsave(filename = "s5s.tiff", plot = s5s, width = 3.5, height = 3, units = "in")


#####



sce_score <- AddModuleScore(
  object = obj_NK,
  features = list(cytotoxic2),
  name = 'cytotoxic_score'
)


metadata2<- sce_score@meta.data[,c("Race_Type","NK_subtype","cytotoxic_score1")]

s5t<- ggplot(metadata2, aes(x = cytotoxic_score1, colour = NK_subtype)) +scale_color_manual(values = color.liberal)+
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

ggsave(filename = "s5t.tiff", plot = s5t, width = 5.4, height = 3, units = "in")



