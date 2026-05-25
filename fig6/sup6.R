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

setwd("~/diskD/lung/paneltotal/fig6/")

obj_B<-qread("obj_B.qs")
obj_Bmem<-qread("obj_Bmem.qs")
obj_PC<-qread("obj_PC.qs")


color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                  '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#####



s6a <- DimPlot(obj_B, reduction = "umap.endo.rpca", label=T, cols = col_vector,group.by = "seurat_clusters",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")


ggsave(filename = " s6a.png", plot =  s6a, width = 4.5, height = 4, units = "in")


shortlist<- c("MS4A1", "FCER2","AICDA","CD80","XBP1","MZB1","SDC1")


plots <- lapply(1:7, function(i) {
  if (i == 5){
    FeaturePlot(obj_B,features = shortlist[i], reduction = "umap.endo.rpca",cols = c('grey','red'))+labs(y = "UMAP_2", x = "UMAP_1",max.cutoff = 4.6)+theme(legend.position = "none")
  }else{
    FeaturePlot(obj_B,features = shortlist[i], reduction = "umap.endo.rpca",cols = c('grey','red'))+labs(y = " ", x = " ",max.cutoff = 4.6,no.axes = TRUE)+theme(legend.position = "none")
  }
})
 s6b <- wrap_plots(plots, ncol = 4) +  plot_layout(guides = 'collect') & 
  guides(fill = guide_legend(title = "Expression"))

ggsave(filename = " s6b.png", plot =  s6b, width = 7.6, height = 3.6,units = "in")



#####

cluster.averages <- AverageExpression(obj_B, return.seurat = T, group.by = 'RNA_snn_res.0.8')

my_order <- c( "g5", 
               "g13", "g14","g17",
               "g0", "g1", "g2", "g9",
               "g12", "g15", "g8", "g16", "g18",
               "g3", "g4", "g6", "g7", "g10", "g11")

B_marker <- c("MS4A1","CD19","BANK1","SELL")

naive_B_markers <- c("TCL1A","FCER2","IL4R")

activated_GC_markers <- c("AICDA","RGS13","CD83","CD86","MYC")

memory_B_markers <- c("AIM2","TNFRSF13B","GPR183","CD80")

plasmablast_markers <- c("MZB1","XBP1","CD38","CD27")

plasma_cell_markers <- c("SDC1","JCHAIN","MZB1","XBP1","CD38","TNFRSF17")

ig_genes <- c("IGHM","IGHD","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHA2")

mhc_ii_genes <- c("HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5")


all_B_markers <- c(B_marker,naive_B_markers,activated_GC_markers,
                   memory_B_markers,plasmablast_markers,plasma_cell_markers,
                   ig_genes,mhc_ii_genes)
all_B_markers <- unique(all_B_markers)

Idents(cluster.averages) <- factor(Idents(cluster.averages), levels = my_order)

 s6c <- DoHeatmap(
  object = cluster.averages,
  features = all_B_markers,  # Genes
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
 s6c
ggsave(filename = " s6c.png", plot =  s6c, width = 6, height = 8, units = "in")

#####



 s6d <- DimPlot(obj_B, reduction = "umap.endo.rpca", label=T, cols = color.liberal,group.by = "stage1_res.0.8",label.size = 4,repel = FALSE) + NoLegend()+
  labs(y = "UMAP_2", x = "UMAP_1")+ggtitle("")


ggsave(filename = " s6d.png", plot =  s6d, width = 4.5, height = 4, units = "in")


#####


pt <- table(obj_B@meta.data$Sample_ID, obj_B@meta.data$stage1_res.0.8)
pt <- as.data.frame(pt)





 s6e1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x =element_blank())+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




 s6e2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))







combined_plot <- ( s6e1 /  s6e2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")

ggsave(filename = " s6e.png", plot = combined_plot, width = 11, height = 5, units = "in")


#####



pt <- table(obj_B@meta.data$T_N, obj_B@meta.data$stage1_res.0.8)
pt <- as.data.frame(pt)




 s6f1 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.9) +
  xlab("") +
  ylab("Fraction")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))




 s6f2 <-  ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "stack", width = 0.9) +
  xlab("") +
  ylab("Counts")  +
  scale_fill_manual(values = color.liberal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(axis.line = element_line(colour = "black"))+
  theme( panel.background = element_blank(),text = element_text(size=12))



combined_plot <- ( s6f1 |  s6f2) + 
  plot_layout(guides = "collect") + 
  plot_annotation() + 
  theme(legend.position = "right")


ggsave(filename = " s6f.png", plot = combined_plot, width = 4.3, height = 3, units = "in")


#####

Idents(obj_PC) <- "Race_TNM"

 s6g <- DEenrichRPlot(
  obj_PC,
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


 s6g[[1]] <-  s6g[[1]] + 
  labs(title =paste0("Black Tumor ","Plasma"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.75) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )


 s6g[[2]] <-  s6g[[2]] + 
  labs(title =paste0("White Tumor ","Plasma"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.75) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )

ggsave(filename = " s6g.png", plot =  s6g, width = 6.5, height = 3, units = "in")


#####


Idents(obj_Bmem) <- "Race_TNM"

 s6h <- DEenrichRPlot(
  obj_Bmem,
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


 s6h[[1]] <-  s6h[[1]] + 
  labs(title =paste0("Black Tumor ","Memory B"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.75) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )


 s6h[[2]] <-  s6h[[2]] + 
  labs(title =paste0("White Tumor ","Memory B"))+
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "yellow",linewidth = 0.75) +theme(
    plot.title  = element_text(size = 14) # Resize title while you're at it
  )

ggsave(filename = " s6h.png", plot =  s6h, width = 6.5, height = 3, units = "in")

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



obj_nB <- subset(obj_B, stage1_res.0.8 == "Naive B")


Idents(obj_nB) <- "Race_TNM"



# Here we run standard DE to get adjusted p-values properly.
markers <- FindMarkers(
  obj_nB,
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
s6i <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw() +
  labs(
    title = "Black Tumor vs White Tumor Naive B",
    x = "",
    y = "Average Log2FC",
    size = "-log10(p-adj)",
    color = "Log2FC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.box = "horizontal")+
  guides( color = "none")

ggsave(filename = "s6i.png", plot = s6i, width = 4.3, height = 3.5, units = "in")


#####



obj_Bact <- subset(obj_B, stage1_res.0.8 == "Activated B")

Idents(obj_Bact) <- "Race_TNM"

# Here we run standard DE to get adjusted p-values properly.
markers <- FindMarkers(
  obj_Bact,
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
s6j <- ggplot(subset_degs, aes(x = gene, y = avg_log2FC)) +
  geom_point(aes(size = -log10(p_val_adj), color = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(range = c(2, 8)) + 
  theme_bw() +
  labs(
    title = "Black Tumor vs White Tumor Memory B",
    x = "",
    y = "Average Log2FC",
    size = "-log10(p-adj)",
    color = "Log2FC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.box = "horizontal")+
  guides( color = "none")

ggsave(filename = "s6j.png", plot = s6j, width = 4.3, height = 3.5, units = "in")


