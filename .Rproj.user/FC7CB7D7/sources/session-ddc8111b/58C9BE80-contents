### Contrastive PCA function (no automatic calculation for alpha)


cPCA <- function(X_target, Y_background, alpha = 1, scale = TRUE, center = TRUE, ncomp = NULL) {
  # Preprocess
  X <- scale(X_target[,-1], center = center, scale = scale)
  Y <- scale(Y_background[,-1], center = center, scale = scale)
  
  # Covariances
  cov_X <- cov(X)
  cov_Y <- cov(Y)
  
  # Contrastive covariance matrix
  cov_c <- cov_X - alpha * cov_Y
  
  # Eigen decomposition
  eig <- eigen(cov_c)
  V <- eig$vectors  # eigenvectors
  
  # Optional limit on components
  if (!is.null(ncomp)) {
    V <- V[, 1:ncomp, drop = FALSE]
  }
  
  # Project datasets
  scores_X <- X %*% V
  scores_Y <- Y %*% V
  
  # Compute per-component variance breakdown
  var_target     <- diag(t(V) %*% cov_X %*% V)
  var_background <- diag(t(V) %*% cov_Y %*% V)
  contrastive_variance <- var_target - alpha * var_background
  
  
  # Return
  list(
    eigenvectors = V,
    eigenvalues = eig$values[1:ncomp], 
    var_target = var_target,
    var_background = var_background,
    contrastive_variance = contrastive_variance,
    scores_target = scores_X,
    scores_background = scores_Y,
    alpha = alpha,
    X_target = X_target,
    Y_background = Y_background,
    center = attr(X, "scaled:center"),
    scale = attr(X, "scaled:scale")
  )
  
  
}