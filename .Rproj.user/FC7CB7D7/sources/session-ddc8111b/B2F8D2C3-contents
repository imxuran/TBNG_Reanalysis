### Heatmap function called in WGCNA function


# Create a color function based on standardized scale
color_func <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#67a9cf", "#f7f7f7", "#ef8a62")
)


make_module_heatmap <- function(module_name, expression_mat, metadata_df,
                                module_gene, module_eigengenes_df) {
  
  
  # Extract a vector of variables IDs that correspond to this module and group
  df_comb <- metadata_df[, c('group', module_gene)]
  
  # Combine the eigengene data and the original data
  df_comb <- cbind(df_comb, module_eigengenes_df[paste0('ME',module_name)])
  
  # Sort by group
  df_comb <- dplyr::arrange(df_comb, group)
  
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- HeatmapAnnotation(
    group = df_comb$group,  # Supply treatment labels
    module_eigengene = anno_barplot(df_comb[[paste0('ME',module_name)]]),  # Add annotation barplot
    col = list(group = c("case" = "#f1a340", "control" = "#998ec3"))  # Pick colors for each experimental group in time_point
  )
  
  # Normalize the original values and keep for heatmap visualization
  mod_mat <- t(scale(df_comb[, -c(1, ncol(df_comb))]))
  
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = paste0('ME',module_name),
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = TRUE,
                                     show_column_names = FALSE
  )
  
  # Return heatmap
  return(heatmap)
}