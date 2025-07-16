### Plot the strongest effect sizes in a forest plot


# Group the families, determine families in each picture
split_list <- list(
  c("Complexity"),
  c("Fucosylation"),
  c(  "Galactosylation", "Sialylation"),
  c("Bisection","Sialylation_L", "Sialylation_E")
)

# List to store the forest plots
forest_list <- list()

# and now follows a buildup of a plot 
for (n in split_list) {
  p <- subset_cohen %>% 
    filter(trait %in% union(pos3,neg3)) %>% 
    filter(family %in% n) %>% 
    ggplot(aes(x=cohen_d, y=trait)) + 
    geom_effect(ggplot2::aes(xmin = CI_lower, xmax = CI_upper, colour = disease, shape = disease), fatten = 6, position = ggstance::position_dodgev(height = 0.5)) + 
    scale_color_manual(values = disease_colors) + 
    scale_shape_manual(values = c(RA = 19, AIH = 19, UC = 17, CD = 15, LMS = 17, COVID19 = 15, MASH = 19, T2D = 17, BC = 15, CRC = 19, MM = 17, PC = 15, PDAC = 17)) +
    theme_forest() + 
    geom_stripes(odd = "#33333333", even = "#00000000") +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", cex = 1, alpha = 0.5) +  
    ggforce::facet_col(facets = ~ family, scales = "free_y", space = "free") + 
    xlim(-4, 4) + xlab("Cohen's D (CI)") + ylab("Glyco traits") 
  
  # As the legends for each subplot might be not identical, here only keep the complete one for all
  if (length(unique((subset_cohen %>% filter(family %in% n))$disease))==13) {
    p <- p + theme( legend.position = "bottom") + guides(color = guide_legend(nrow = 1, override.aes = list(size = 3)))
  } else {
    p <- p + guides(color = "none", fill = "none", shape = "none")
  }
  
  forest_list <- append(forest_list,list(p))
}

## Merge plots

# Merge 4 forest plots stored in the list with patchwork
p <- (wrap_plots(forest_list, ncol=4) + plot_annotation(tag_levels = 'A')) + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") 
#& guides(color = guide_legend(nrow = 1, override.aes = list(size = 3)), shape = guide_legend(nrow = 1, override.aes = list(size = 3)))

# Save the combined plot
ggsave("Results/EffectSize_ForestPlot.png", plot = p, width = 12, height = 8 , units = "in", dpi = 500)
