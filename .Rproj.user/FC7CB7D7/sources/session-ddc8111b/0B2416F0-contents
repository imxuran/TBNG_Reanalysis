### function for plotting scores and loadings


cPCA_plot <- function(alpha){
  
  # Conduct cPCA
  res <- cPCA(X_target = df_all_case, Y_background = df_all_control, ncomp=10, alpha = alpha)
  
  # Plot the score figure
  # Only select the top five PCs
  df_plot <- rbind(
    data.frame(res$scores_target[, 1:5], group = "case") %>% bind_cols(df_all_case %>% select(disease)),
    data.frame(res$scores_background[, 1:5], group = "control")%>% bind_cols(df_all_control %>% select(disease))
  )
  
  names(df_plot)[1:5] <- c("cPC1", "cPC2", "cPC3", "cPC4", "cPC5")
  
  levels(df_plot$disease) <- c(disease_order,"control")
  
  df_plot[df_plot$group=='control',]$disease <- 'control'
  
  df_plot$disease <- factor(df_plot$disease, levels = c(disease_order,'control'))
  
  # extract loadings
  loadings <- as.data.frame(res$eigenvectors[, 1:2])  # First 2 PCs for variables
  colnames(loadings) <- c('cPC1','cPC2')
  rownames(loadings) <- colnames(df_all_case[,-1])
  loadings$Variable <- rownames(loadings)
  loadings$Family <- NA
  for (family in names(families)) {
    loadings$Family[loadings$Variable %in% families[[family]]] <- family
  }
  loadings$length <- sqrt(loadings[,1]^2 + loadings[,2]^2) # Add loadings' magnitude

  
  
  # Score plot
  p <- ggplot(df_plot, aes(x = cPC1, y = cPC2, color = disease)) +
    geom_point(size = 1) +
    labs(title = paste("cPCA (alpha=",alpha,")"), x = "cPC1", y = "cPC2") +
    scale_color_manual(values=c(disease_colors, control = '#FEF295'))+
    theme_minimal() + stat_ellipse() + theme( legend.position = "bottom") + guides(color = guide_legend(nrow = 1))
  
  p <- ggMarginal(p, type = 'density', groupColour = T, groupFill = T)
  ggsave(paste0("Results/cPCA/cPCA_scores(alpha=",alpha,").png"), p, width = 14, height = 8, dpi = 300)
  
  
  # Calculate statistics index of scores for later centroids visualization
  summary_stats <- df_plot %>%
    group_by(disease) %>%
    summarise(x_mean = mean(cPC1), y_mean = mean(cPC2), 
              x_min = min(cPC1), x_max = max(cPC1),
              y_min = min(cPC2), y_max = max(cPC2),
              x_ci_lower = mean(cPC1) - 1.96 * sd(cPC1) / sqrt(n()),
              x_ci_upper = mean(cPC1) + 1.96 * sd(cPC1) / sqrt(n()),
              y_ci_lower = mean(cPC2) - 1.96 * sd(cPC2) / sqrt(n()),
              y_ci_upper = mean(cPC2) + 1.96 * sd(cPC2) / sqrt(n()),
              x_se = sd(cPC1) / sqrt(n()), 
              y_se = sd(cPC2) / sqrt(n()), .groups = 'drop')
  
  
  # Scores with centroids
  p <- ggplot() +
    
    geom_point(data = df_plot, aes(x = cPC1, y = cPC2, color = disease), 
               alpha = 0.1, size = 1.5) +
    
    geom_point(data = summary_stats, aes(x = x_mean, y = y_mean, color = disease), 
               size = 4, shape = 19, show.legend = FALSE) +
    
    geom_text_repel(data = summary_stats,
                    aes(x = x_mean, y = y_mean, label = disease, color=disease), 
                    #color = loadings$Arrow_Color,
                    size = 4,
                    alpha = 1,
                    max.overlaps = 100, inherit.aes = FALSE,
                    show.legend = FALSE) + 
    
    geom_errorbar(data = summary_stats, 
                  aes(x = x_mean, y = y_mean, 
                      xmin = x_mean - 10*x_se, xmax = x_mean + 10*x_se, color = disease),
                  width = 0.1) +
    
    geom_errorbar(data = summary_stats, 
                  aes(x = x_mean, y = y_mean, 
                      ymin = y_mean - 10*y_se, ymax = y_mean + 10*y_se, color = disease),
                  width = 0.1) +
    
    labs(title = paste("cPCA (alpha=",alpha,")"), x = "cPC1", y = "cPC2") +
    scale_color_manual(values=c(disease_colors, control = '#FEF295'))+
    theme_minimal() + stat_ellipse() + theme( legend.position = "bottom") + guides(color = guide_legend(nrow = 1))
  
  p <- ggMarginal(p,type = 'density', groupColour = T, groupFill = T)
  
  ggsave(paste0("Results/cPCA/scores(alpha=",alpha,").png"), p, width = 14, height = 8, dpi = 300)
  
  # Loadings magnitude lolliplot
  p <- ggplot(loadings, aes(x = reorder(Variable, length), y = length, color = Family)) +
    geom_segment(aes(xend = Variable, y = 0, yend = length), linewidth = 0.7) +
    geom_point(size = 3) +
    coord_flip() +
    theme_minimal() +
    labs(x = NULL, y = NULL , title = "Loadings magnitude", color = "Family") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      panel.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(), 
      axis.text.x = element_blank(),  
      axis.ticks.x = element_blank(),  
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_markdown(size = 10)
    ) + scale_color_manual(values=family_colors)
  
  ggsave(paste0("Results/cPCA/loading_magnitude(alpha=",alpha,").png"), p, width = 5, height = 6, dpi = 300, bg = "white")
  
}