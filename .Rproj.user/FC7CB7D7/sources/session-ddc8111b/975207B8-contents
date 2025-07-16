plot_list <- list(
  Fucosylation = c('A3F','A3Fa','CFa'),
  Sialylation = c('A3E','A4E','A3L','A4L'),
  Bi_antennaty = c('A2G','A2E','A2L','A2F','CB'),
  Complexity = c('THy','TM','TC','MM')
)


for (n in names(plot_list)) {
  
  data <- dmfa_loadings[dmfa_loadings$Variable %in% plot_list[[n]],]
  
  p <- ggplot(data, aes(x = Dim1, y = Dim2, color = Group)) +
    
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    
    geom_point(data = data, size = 1.5, aes(x = Dim1, y = Dim2,color= Group,shape = Variable), alpha = 1) + 
    
    geom_text_repel(data = data,
                    aes(x = Dim1, y = Dim2, label = Variable, color=Group),
                    size = 3,
                    alpha = 1,
                    max.overlaps = 100,
                    show.legend = FALSE) +
    
    theme_minimal() +
    theme(legend.position = "right") +
    labs(
      title = paste(n,' Traits'),
      x = sprintf("Dimension %d (%.1f%%)", 1, prop_var[1]),
      y = sprintf("Dimension %d (%.1f%%)", 2, prop_var[2]),
      color = 'group'
    ) + scale_color_manual(values=disease_colors) + scale_shape_discrete()

  ggsave(paste0("Results/DMFA/DMFA_Loadings_",n,".png"), p, width = 7, height = 5, dpi = 300, bg = "white")
  
}

