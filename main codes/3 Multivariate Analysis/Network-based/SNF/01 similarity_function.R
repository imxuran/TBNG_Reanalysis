### Compute similarity matrix


# Functions to compute similarity
compute_similarity <- function(data, method = "correlation", non_negative = TRUE) {
  if (method == "correlation") {
    corr_matrix <- cor(data, use = "pairwise.complete.obs")
    
    if (non_negative) {
      similarity_matrix <- (corr_matrix + 1) / 2  # Scale to [0, 1]
    } else {
      similarity_matrix <- corr_matrix 
    }
    #similarity_matrix <- 1 - similarity_matrix
    
  } else if (method == "euclidean") {
    
    dist_matrix <- as.matrix(dist(t(data), method = "euclidean"))  
    similarity_matrix <- 1 / (1 + dist_matrix) 
    
  } else if (method == "mahalanobis") {
    
    dist_matrix <- as.matrix(vegdist(t(data), method = "mahalanobis"))
    
    similarity_matrix <- 1 / (1 + dist_matrix)
    
    
  } else {
    stop("Unknown method. Choose 'correlation' or 'euclidean'")
  }
  
  
  # Reorder the columns and rows of the matrix based on the new order defined above
  similarity_matrix <- similarity_matrix[intersect(trait_order,rownames(similarity_matrix)),intersect(trait_order,rownames(similarity_matrix))]
  
  return(similarity_matrix)
}