library(ggplot2)
library(dplyr)

library(monocle)
library(igraph)
library(cowplot)
library(qs)



color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0',
                  '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8',
                  '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')


setwd("~/scRNA/myeloid")

myeloid_traj <- qread("myeloid_traj.qs")


#####

p6 <- plot_cell_trajectory(myeloid_traj, color_by = "State", cell_size = 0.02, show_branch_points = T) + 
  scale_color_manual(values = color.liberal)

p3 <- plot_cell_trajectory(myeloid_traj, color_by = "T_N", show_branch_points = F, cell_size = 0.05) +
  scale_color_manual(values = color.liberal) 

p4 <- plot_cell_trajectory(myeloid_traj, color_by = "cell_type", show_branch_points = F, cell_size = 0.05) +
  scale_color_manual(values = color.liberal)


p5 <- plot_cell_trajectory(myeloid_traj, color_by = "R_T_N", show_branch_points = F, cell_size = 0.05) + 
  scale_color_manual(values = color.liberal)



f4g <- plot_grid(p6+scale_x_reverse(), p3+scale_x_reverse(), p4+scale_x_reverse(), p5+scale_x_reverse(),  ncol = 4)


ggsave(filename = "f4g.png", plot = f4g, width = 15, height = 4, units = "in")



#####


s4k <- plot_cell_trajectory(myeloid_traj, color_by = "cell_type", cell_size = 0.03) + 
  scale_color_manual(values = color.liberal) +
  facet_grid(R_T_N ~ cell_type)+
  theme(
    legend.position = "none",
  )

ggsave(filename = "s4k.png", plot = s4k, width = 7, height = 8, units = "in")


