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


color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


setwd("~/diskD/lung/paneltotal/fig1/")

rds1 <-qread("rds1.qs")


#####


s1b <- DimPlot(rds1, reduction = "umap.rpca", label=T, cols = col_vector,group.by = "seurat_clusters",label.size = 4) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")
s1b
ggsave(filename = "s1b.png", plot = s1b, width = 5, height = 4.5, units = "in")






#####

shortlist <- c("MS4A1","PECAM1","EPCAM","DCN","TPSB2","LYZ","MZB1","KLRF1","CD3D")

plots <- lapply(1:9, function(i) {
  if (i == 9){
    FeaturePlot(rds1,features = shortlist[i], reduction = "umap.rpca",cols = c('grey','red'))+labs(y = "UMAP_2", x = "UMAP_1",max.cutoff = 4.6)+theme(legend.position = "none")
  }else{
    FeaturePlot(rds1,features = shortlist[i], reduction = "umap.rpca",cols = c('grey','red'))+labs(y = " ", x = " ",max.cutoff = 4.6,no.axes = TRUE)+theme(legend.position = "none")
  }
})
s1c <- wrap_plots(plots, ncol = 2) +  plot_layout(guides = 'collect') & 
  guides(fill = guide_legend(title = "Expression"))
#  theme(legend.position = "right") &

ggsave(filename = "s1c.png", plot = s1c, width = 3.9, height = 9,units = "in")


#####



s1d <- DimPlot(rds1, reduction = "umap.rpca", label=T, cols = color.liberal,group.by = "celltype2_RNA_snn_res.2",label.size = 5,repel = TRUE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")

ggsave(filename = "s1d.png", plot = s1d, width = 5, height = 4.5, units = "in")


#####

b_cell_markers <- c("CD79A", "CD79B", "CD22", "CD19", "MS4A1", "PAX5", "BANK1")

endothelial_markers <- c("CAVIN2", "GNG11", "CLDN5", "VWF", "RAMP2", "ENG", 
                         "ADGRL4", "PECAM1", "EMCN")

epithelial_markers <- c("AKR1C1", "KRT8", "KRT18", "KRT19", "EPCAM", "CD24", 
                        "MUC1", "CDH1")

fibroblast_markers <- c("IGFBP7", "CALD1", "FSTL1", "CCN2", "THY1", "COL6A3", 
                        "LUM", "DCN", "COL1A2")

mast_cell_markers <- c("FCER1A", "TPSAB1", "MS4A2", "RGS13", "HPGDS", "HDC", 
                       "IL1RL1", "KIT")

myeloid_markers <- c("S100A8", "S100A9", "LYZ", "CD74", "HLA-DRA", "IL1B", 
                     "CD14", "CD33", "ITGAM", "CSF1R")

nk_markers <- c("GZMH", "KLRB1", "CST7", "CTSW", "CCL5", "NKG7", "GNLY", 
                "KLRD1", "SPON2", "FGFBP2")

plasma_markers <- c("IGKC", "IGLC2", "IGHG1", "IGLC3", "JCHAIN", "IGHG2", "MZB1")

t_cell_markers <- c("CD3D", "CD3E", "CD27", "BATF", "IL7R", "CD2", "CD8A", 
                    "CD3G", "TRAC")
all_markers <- c(b_cell_markers,endothelial_markers,epithelial_markers,fibroblast_markers,
                 mast_cell_markers,myeloid_markers,nk_markers,plasma_markers,t_cell_markers)





b_cell_clusters <- paste0("g", c(1, 52))
endothelial_clusters <- paste0("g", c(47, 33))
epithelial_clusters <- paste0("g", c(39, 28, 13, 34, 48, 20, 37, 42, 32, 40))
fibroblast_clusters <- paste0("g", c(17, 41, 44))
mast_cell_clusters <- paste0("g", c(54, 29, 53, 50))
myeloid_clusters <- paste0("g", c(43, 16, 27, 26, 49, 10, 36, 14, 12, 9, 35, 19, 18, 21))
nk_clusters <- paste0("g", c(7, 5, 8, 6))
plasma_clusters <- paste0("g", c(46, 25, 31, 45, 51))
t_cell_clusters <- paste0("g", c(30,55, 56, 4, 3, 2, 24, 22, 23, 15, 11, 0, 38, 57))

ordered_clusters <- c(
  b_cell_clusters, endothelial_clusters, epithelial_clusters, 
  fibroblast_clusters, mast_cell_clusters, myeloid_clusters, 
  nk_clusters, plasma_clusters, t_cell_clusters
)


cluster.averages <- AggregateExpression(rds1, return.seurat = T, group.by = 'seurat_clusters')
Idents(cluster.averages) <- factor(Idents(cluster.averages) , levels =ordered_clusters )


s1e <- DoHeatmap(
  object = cluster.averages,
  features = all_markers, 
  size = 5,
  draw.lines = FALSE,
  disp.min = 0,
  disp.max = 2,
) +
  scale_fill_gradientn(colors = c("lightgrey", "blue")) +
  theme(
    axis.text.y = element_text(size = 10, family = "Arial", color = "black"),  # Gene labels
    axis.text.x.bottom = element_text(size = 8, angle = 90,hjust = 1, family = "Arial", color = "black"),  # Bottom axis cluster labels
    axis.text.x.top = element_blank(),  #axis.text.x.top = element_text(size = 8, angle = 45,  family = "Arial",hjust = 1, family = "Arial", color = "black"),  # Top axis cluster labels
    legend.position = "none",  # Removes the legend
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank()   # Remove y-axis title
  )
ggsave(filename = "s1e.png", plot = s1e, width = 9, height = 12, units = "in")



#####



pt <- table(rds1@meta.data$Sample_ID, rds1@meta.data$celltype2_RNA_snn_res.2)
pt <- as.data.frame(pt)


s1f1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x =element_blank())+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=15))




s1f2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=15))







combined_plot <- (s1f1 / s1f2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")
combined_plot

ggsave(filename = "s1f.png", plot = combined_plot, width = 14, height = 7.5, units = "in")
