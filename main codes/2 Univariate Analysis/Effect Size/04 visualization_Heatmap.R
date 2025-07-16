### Plot all effect sizes in a heatmap along with hierarchical clustering


# Bind all effect sizes together
cohen_df <- do.call(cbind,lapply(cohen_list,function(x) x[intersect(trait_order, union(pos3,neg3)),2]))

rownames(cohen_df) <- intersect(trait_order, union(pos3,neg3))

# Considering there are traits absent in some sets, set the null value to zero.
cohen_df[is.na(cohen_df)] <- 0

# Add annotation of the trait and disease group for later plot
ann_row_group <- sapply(rownames(cohen_df), function(x) {
  names(families)[which(sapply(families, function(f) x %in% f))]
})
ann_row <- as.data.frame(ann_row_group)
colnames(ann_row) <- c("Trait")

ann_col_group <- sapply(colnames(cohen_df), function(x) {
  names(disease_types)[which(sapply(disease_types, function(f) x %in% f))]
})
ann_col <- as.data.frame(ann_col_group)
colnames(ann_col) <- c("Disease")


ann_colors = list(
  Trait = family_colors,
  Disease = type_colors
)

# Optional: set zero as the median color
data_range <- max(abs(c(min(cohen_df, na.rm = TRUE), 
                        max(cohen_df, na.rm = TRUE))))
breaks <- seq(-data_range, data_range, length.out = 101)


png("Results/EffectSize_Heatmap.png", width = 3200, height = 2400, res = 300)

pheatmap::pheatmap(cohen_df, cluster_rows = TRUE, cluster_cols = TRUE,
                   show_rownames = TRUE, show_colnames = TRUE, 
                   annotation_row = ann_row,
                   annotation_col = ann_col,
                   annotation_colors = ann_colors, 
                   border_color = NA,
                   #breaks = breaks,
                   clustering_distance_rows = "correlation",
                   main = 'Effect Size Heatmap')

dev.off()