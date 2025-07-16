### Function of consensus WGCNA (to detect the shared modules in multiple datasets)

consensus_wgcna <- function(df_list, n) {
  
  df_raw_list <- df_list
  
  for (i in 1:length(df_list)) {
    df_list[[i]] <- list(
      data = as.data.frame(scale(df_list[[i]][, analy_traits]))
    )
  }
  
  consensusModules <- blockwiseConsensusModules(df_list, power = 6, minModuleSize = 5,  networkType = "unsigned")
  
  png(paste0("Results/WGCNA/consensus/",n,"_dendrogram.png"), width = 2500, height = 1600, res = 300)
  
  plotDendroAndColors(
    dendro = as.hclust(consensusModules$dendrograms[[1]]),
    colors = cbind(consensusModules$colors),
    groupLabels = "Consensus modules",
    main = paste("Consensus gene dendrogram and module colors -", n),
    dendroLabels = analy_traits,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  
  dev.off()
  
  consensus_kME <- consensusKME(df_list, moduleLabels = consensusModules$colors, multiEigengenes = consensusModules$multiMEs)
  
  
  ### Which modules have biggest differences across treatment groups
  for (n in names(df_raw_list)) {
    
    df_raw <- df_raw_list[[n]]
    module_eigengenes <- consensusModules$multiMEs[[n]]$data
    dynamicColors <- consensusModules$colors
    
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
      module_genes <- names(dynamicColors)[dynamicColors == substr(top_module, 3, nchar(top_module))]
      
      
      # a visual support
      p1 <- ggplot(ME_df, aes(x = group, y = !!sym(top_module), color = group)) + geom_boxplot(width = 0.2, outlier.shape = NA) + ggforce::geom_sina(maxwidth = 0.3) + theme_classic()
      
      p2 <- make_module_heatmap(module_name = substr(top_module, 3, nchar(top_module)), expression_mat = df, metadata_df = df_raw, module_gene = module_genes, module_eigengenes_df = module_eigengenes)
      
      png(paste0("Results/WGCNA/consensus/",n,"_heatmap.png"), width = 4500, height = 2300, res = 300)
      
      grid.arrange(p1, grid.grabExpr(draw(p2)), ncol = 2,
                   top = textGrob(n, gp = gpar(fontsize = 16, font = 2)))
      
      dev.off()
      
    }
    
  }
  
  return(consensus_kME)
}


## cWGCNA implementation
conw_C <- consensus_wgcna(data_list[names(data_list) %in% disease_types$C], 'C')

conw_AI <- consensus_wgcna(data_list[names(data_list) %in% disease_types$AI], 'AI')

conw_INF <- consensus_wgcna(data_list[names(data_list) %in% disease_types$INF], 'INF')

conw_M <- consensus_wgcna(data_list[names(data_list) %in% disease_types$M], 'M')

consensus_wgcna(data_list, 'All')

consensus_wgcna(data_list[names(data_list) %in% c('CRC','PC','BC','PDAC')], 'Solid C')