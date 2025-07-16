### Function for WGCNA


wgcna <- function(df, n) {
  
  df_raw <- df

  df <- scale(df[, analy_traits])
  
  #rownames(df) <- df_raw$SampleID
  
  ### Determine Soft Thresholding Power
  # Choosing a soft-threshold to fit a scale-free topology to the network
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  
  # Call the network topology analysis function
  sft=pickSoftThreshold(df, dataIsExpr = TRUE, powerVector = powers, corFnc = bicor, networkType = "unsigned")
  
  softPower = which.max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
  
  
  # Plot the soft cutoff results
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",type="n", 
       main = paste("Scale independence -",n))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       labels=powers, cex=0.9, col="red")
  
  abline(h=0.80,col="red")  # Red line corresponds to using an R^2 cut-off
  
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity -",n))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
  
  
  
  ### Generating adjacency and TOM similarity matrices based on the selected softpower
  #calclute the adjacency matrix
  adj= adjacency(df, type = "unsigned", power = softPower)
  
  #turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
  TOM=TOMsimilarity(adj)
  
  colnames(TOM) = rownames(TOM) = colnames(df)
  dissTOM=1-TOM
  
  #hierarchical clustering of the genes based on the TOM dissimilarity measure
  geneTree = hclust(as.dist(dissTOM), method = "average")
  
  # Module identification using dynamic tree cut
  dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = 5)
  dynamicColors = labels2colors(dynamicMods)
  
  png(paste0("Results/WGCNA/individual/",n,"_dendrogram.png"), width = 2500, height = 1600, res = 300)
  
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = NULL, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = paste("Dendrogram and module colors -",n))
  
  dev.off()
  
  
  
  ### Which modules have biggest differences across treatment groups
  module_eigengenes <- moduleEigengenes(df, colors = dynamicColors)$eigengenes
  # Create the design matrix from the `group` variable
  des_mat <- model.matrix(~ df_raw$group)
  fit <- limma::lmFit(t(module_eigengenes), design = des_mat)
  # Apply empirical Bayes to smooth standard errors, very practical in our case
  fit <- limma::eBayes(fit)
  # Apply multiple testing correction and obtain stats in a single object
  stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>% tibble::rownames_to_column("module")
  
  # exploring moduls
  ME_df <- cbind(module_eigengenes,df_raw['group'])
  
  stats_df <- stats_df %>%
    dplyr::filter(module != "ME0" & module != "MEgrey") %>%
    dplyr::arrange(adj.P.Val)
  
  if (nrow(stats_df) == 0) {
    # Case: Only ME0 exists
    print("No modules identified!")
  } else {
    top_module <- stats_df$module[1]
    
    # Get a vector of variables IDs that correspond to this (most significant) module
    module_genes <- colnames(df)[dynamicColors == substr(top_module, 3, nchar(top_module))]
    
    
    # a visual support
    p1 <- ggplot(ME_df, aes(x = group, y = !!sym(top_module), color = group)) + geom_boxplot(width = 0.2, outlier.shape = NA) + ggforce::geom_sina(maxwidth = 0.3) + theme_classic()
    
    p2 <- make_module_heatmap(module_name = substr(top_module, 3, nchar(top_module)), expression_mat = df, metadata_df = df_raw, module_gene = module_genes, module_eigengenes_df = module_eigengenes)
    
    png(paste0("Results/WGCNA/individual/",n,"_heatmap.png"), width = 4500, height = 2300, res = 300)
    
    grid.arrange(p1, grid.grabExpr(draw(p2)), ncol = 2,
                 top = textGrob(n, gp = gpar(fontsize = 16, font = 2)))
    
    dev.off()
    
    return(module_genes)
  }
}
