### SNF computation function


# Functions to construct the consensus similarity matrix
compute_SNF <- function(similarity_list) {
  K = 7		# number of neighbors
  alpha = 0.5  	# hyperparameter
  T = 10 	# Number of Iterations
  
  #truelabel = list() ##the ground truth of the simulated data;
  
  #similarity_list = lapply(similarity_list, standardNormalization)
  
  wlist <- list()  # Store similarity graphs
  
  for (d in similarity_list) {
    #d <- dist2(as.matrix(t(s)),as.matrix(t(s)))
    w <- affinityMatrix(d, K, alpha) # Computes the affinity matrix for a given distance matrix
    wlist <- append(wlist, list(w))
  }
  
  
  W <- SNF(wlist, K, T)  # Compute the overall network
  
  C = 5		# number of clusters
  group <- spectralClustering(W, C) 	# the final subtypes information
  
  col = getColorsForGroups(group, colors = colorRampPalette(brewer.pal(9, "Paired"))(9))
  
  displayClustersWithHeatmap(W, group, col = colorRampPalette(brewer.pal(9, "Purples"))(50), ColSideColors=col)
  
  return(W)
}


### Construct the consensus similarity matrix

## For case samples

# total consensus matrix
png("Results/SNF/case_total_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_total_case <- compute_SNF(similarity_case)
dev.off()

# consensus matrix of cancer
png("Results/SNF/case_cancer_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_cancer_case <- compute_SNF(similarity_case[names(similarity_case) %in% disease_types$C])
dev.off()

# consensus matrix of autoimmune & inflammatory
png("Results/SNF/case_ai_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_ai_case <- compute_SNF(similarity_case[names(similarity_case) %in% disease_types$AI])
dev.off()

# consensus matrix of infectious
png("Results/SNF/case_inf_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_inf_case <- compute_SNF(similarity_case[names(similarity_case) %in% disease_types$INF])
dev.off()

# consensus matrix of metabolic
png("Results/SNF/case_meta_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_meta_case <- compute_SNF(similarity_case[names(similarity_case) %in% disease_types$M])
dev.off()


### Construct the consensus similarity matrix

## For control samples
# total consensus matrix for controls
png("Results/SNF/control_total_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_total_control <- compute_SNF(similarity_control)
dev.off()

# consensus matrix of cancer
png("Results/SNF/control_cancer_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_cancer_control <- compute_SNF(similarity_control[names(similarity_control) %in% disease_types$C])
dev.off()

# consensus matrix of autoimmune & inflammatory
png("Results/SNF/control_ai_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_ai_control <- compute_SNF(similarity_control[names(similarity_control) %in% disease_types$AI])
dev.off()

# consensus matrix of infectious
png("Results/SNF/control_inf_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_inf_control <- compute_SNF(similarity_control[names(similarity_control) %in% disease_types$INF])
dev.off()

# consensus matrix of metabolic
png("Results/SNF/control_meta_heatmap.png", width = 2000, height = 2000, res = 500)
consensus_meta_control <- compute_SNF(similarity_control[names(similarity_control) %in% disease_types$M])
dev.off()