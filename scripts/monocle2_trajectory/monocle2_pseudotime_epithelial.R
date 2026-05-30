library(ggplot2)
library(dplyr)

library(monocle)
library(igraph)
library(cowplot)
library(qs)


setwd("~/scRNA/epi")

color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0',
                  '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8',
                  '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

epi_traj.qs <- qread("epi_traj.qs")

#####
p8 <- plot_cell_trajectory(epi_traj.qs, color_by = "State_merged", cell_size = 0.02, show_branch_points = F) + scale_color_manual(values = color.liberal)
p2 <- plot_cell_trajectory(epi_traj.qs, color_by = "Pseudotime", show_branch_points = F, cell_size = 0.05)  + scale_color_viridis_c()
p4 <- plot_cell_trajectory(epi_traj.qs, color_by = "cell_type", show_branch_points = F, cell_size = 0.05) + scale_color_manual(values = color.liberal)
p5 <- plot_cell_trajectory(epi_traj.qs, color_by = "R_T_N", show_branch_points = F, cell_size = 0.05) + scale_color_manual(values = color.liberal)

f2g <- plot_grid(p8, p4, p2, p5,  ncol = 4)
ggsave(filename = "f2g.png", plot = f2g, width = 16, height = 4, units = "in")



#####
s2l <- plot_cell_trajectory(epi_traj.qs, color_by = "cell_type",show_branch_points = F, cell_size = 0.2) + 
  scale_color_manual(values = color.liberal) +
  facet_wrap(~cell_type, nrow = 1)  +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 15)
  )

ggsave(filename = "s2l.png", plot = s2l, width = 12, height = 2, units = "in")
