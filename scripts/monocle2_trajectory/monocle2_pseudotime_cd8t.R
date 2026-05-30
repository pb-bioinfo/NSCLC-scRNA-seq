library(ggplot2)
library(dplyr)

library(monocle)
library(igraph)
library(cowplot)
library(qs)



color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0',
                  '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8',
                  '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')


setwd("~/scRNA/CD8")

CD8_traj <- qread("CD8total.qs")

p6 <- plot_cell_trajectory(CD8_traj, color_by = "State", cell_size = 0.02, show_branch_points = T) +
  scale_color_manual(values = color.liberal)

p2 <- plot_cell_trajectory(CD8_traj, color_by = "Pseudotime", show_branch_points = T, cell_size = 0.05)  +
  scale_color_viridis_c()

p4 <- plot_cell_trajectory(CD8_traj, color_by = "stage2", show_branch_points = F, cell_size = 0.05) + 
  scale_color_manual(values = color.liberal)


p5 <- plot_cell_trajectory(CD8_traj, color_by = "R_T_N", show_branch_points = F, cell_size = 0.05) + 
  scale_color_manual(values = color.liberal)


f3g <- plot_grid(p6, p2, p4, p5,  ncol = 4)

ggsave(filename = "f3g.png", plot = f3g, width = 15, height = 4, units = "in")


#####

s3l <- plot_cell_trajectory(CD8_traj, color_by = "stage2", cell_size = 0.03) + 
  scale_color_manual(values = color.liberal) +
  facet_grid(R_T_N ~ stage2)+
  theme(
    legend.position = "none",
  )

ggsave(filename = "s3l.png", plot = s3l, width = 9, height = 8, units = "in")





