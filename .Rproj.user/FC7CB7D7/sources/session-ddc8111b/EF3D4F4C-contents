### gPCA conduction and results extraction


# Emply PCA
pca_result <- prcomp(data_numeric, scale = TRUE)

# Print PCA summary
# summary(pca_result)

# Extract scores and loadings for visualizaiton
scores <- as.data.frame(pca_result$x) 
scores$Group <- group  # Add group information back
scores$Group <- factor(scores$Group, levels = c(disease_order,'control'))

loadings <- as.data.frame(pca_result$rotation)
loadings$Variable <- rownames(loadings)
loadings$Family <- NA
for (family in names(families)) {
  loadings$Family[loadings$Variable %in% families[[family]]] <- family
}

# scale the loadings length for better visualization
pc_range <- max(abs(range(c(scores$PC1, scores$PC2))))
loading_range <- max(abs(range(c(loadings[,1], loadings[,2]))))
auto_scale <- pc_range / loading_range * 0.3


# Calculate the variances
variances <- (pca_result$sdev^2)/sum(pca_result$sdev^2)

colnames(scores)[1:20] <- paste(colnames(scores[1:20]), sprintf("%.2f%%", variances * 100), sep = ", ")


# Visualize the scores with different PCs
panel_custom <- function(x, y) {
  points(x, y, pch = 20, cex = 0.5, col=c(disease_colors, control = '#FEF295')[as.character(scores[,21])])
}

png("Results/cPCA/gPCA_scores.png", width = 6000, height = 4000, res = 300)
pairs(scores[,1:5], panel = panel_custom,upper.panel = NULL)
dev.off()