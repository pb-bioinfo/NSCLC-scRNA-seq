library(ggplot2)
library(dplyr)
library(patchwork)

library(monocle)
library(igraph)
library(cowplot)
library(qs)

library(tidyverse) 

color.liberal = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0',
                  '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8',
                  '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')


setwd("~/scRNA/Tcell")

CD4_traj <- qread("CD4_traj.qs")

#####

p6 <- plot_cell_trajectory(CD4_traj, color_by = "State", cell_size = 0.02, show_branch_points = T) +
  scale_color_manual(values = color.liberal)

p2 <- plot_cell_trajectory(CD4_traj, color_by = "Pseudotime", show_branch_points = T, cell_size = 0.05)  +
  scale_color_viridis_c()

p4 <- plot_cell_trajectory(CD4_traj, color_by = "stage2", show_branch_points = F, cell_size = 0.05) + 
  scale_color_manual(values = color.liberal)


p5 <- plot_cell_trajectory(CD4_traj, color_by = "R_T_N", show_branch_points = F, cell_size = 0.05) + 
  scale_color_manual(values = color.liberal)


s3p <- plot_grid(p6, p2, p4, p5,  ncol = 4)

ggsave(filename = "f3g.png", plot = s3p, width = 15, height = 4, units = "in")


#####


cell_tensor = table(CD4_traj$stage2,CD4_traj$State,CD4_traj$Race_TNM)
unique_types <- unique(CD4_traj$stage2)



create_cell_proportion_plot <- function(cell_type_name, cell_tensor_data) {
  
  cell_type_slice <- cell_tensor_data[cell_type_name, , c("B_T", "W_T")]
  
  cell_type_df <- as.data.frame.matrix(cell_type_slice) %>%
    rownames_to_column(var = "State") %>%
    pivot_longer(
      cols = c("B_T", "W_T"),
      names_to = "Race_TNM",
      values_to = "Cell_Count"
    )
  cell_proportion_df <- cell_type_df %>%
    group_by(Race_TNM) %>% 
    mutate(Proportion = Cell_Count / sum(Cell_Count)) %>% 
    ungroup() 
  
  p <- ggplot(cell_proportion_df, aes(x = State, y = Proportion, fill = Race_TNM)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("B_T" = "#1f77b4", "W_T" = "#ff7f0e")) +
    labs(
      title = cell_type_name, # Use the cell type name as the title
      x = NULL,
      y = "Proportion of Cells",
      fill = "Patient Group"
    ) +
    theme_classic() +
    theme_bw(base_size = 15) +
    theme(axis.line = element_line(colour = "black")) +
    theme(panel.background = element_blank(), text = element_text(size = 12))
  
  return(p)
}



plot_list2 <- list()

for (type in unique_types) {
  plot_list2[[type]] <- create_cell_proportion_plot(type, cell_tensor)
}

s3r <- wrap_plots(plot_list2, ncol = 4) +
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'right') 


ggsave("s3r.png",plot = s3r, width = 12, height = 4, units = "in")