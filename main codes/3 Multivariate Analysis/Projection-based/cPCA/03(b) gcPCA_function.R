### Generalized contrastive PCA function (automatic calculation for alpha)


# Generalized contrastive PCA (gcPCA) function
gcPCA <- function(Ra, Rb, method = 'v4', Ncalc = Inf, Nshuffle = 0, normalize_flag = TRUE,
                  alpha = 1, alpha_null = 0.975, cond_number = 1e13, verbose = TRUE) {
  
  # Internal function to normalize data using L2 norm (Euclidean)
  normalize <- function(data) {
    data <- scale(data) # Center and scale columns
    data <- data / norm(data, type = "2") # Normalize overall L2 norm
    return(data)
  }
  
  inspect_inputs <- function(Ra, Rb) {
    if (ncol(Ra) != ncol(Rb)) stop("Ra and Rb must have the same number of features")
    
    if (normalize_flag) {
      Ra <- normalize(Ra)
      Rb <- normalize(Rb)
    }
    
    # Determine how many gcPCs can be computed (limited by data rank)
    n_gcpcs <- min(nrow(Ra), nrow(Rb))
    RaRb <- rbind(Ra, Rb)
    svd_result <- svd(RaRb)
    Sab <- svd_result$d
    tol <- max(dim(RaRb)) * .Machine$double.eps * max(Sab)
    
    if (sum(Sab > tol) < n_gcpcs) {
      warning("Input data is rank-deficient! Discarding dimensions.")
      n_gcpcs <- sum(Sab > tol)
    }
    
    # Basis vectors for projection (top right singular vectors)
    J <- svd_result$v[, 1:n_gcpcs, drop = FALSE]
    
    if (method %in% c('v2.1', 'v3.1', 'v4.1')) {
      n_gcpcs <- min(Ncalc, n_gcpcs)
      if (verbose) message(n_gcpcs, " gcPCs will be returned.")
    }
    
    list(n_gcpcs = n_gcpcs, J = J, Ra = Ra, Rb = Rb)
  }
  
  # Core gcPCA computation
  fit <- function(Ra, Rb) {
    inspected <- inspect_inputs(Ra, Rb)
    Ra <- inspected$Ra; Rb <- inspected$Rb
    J <- inspected$J
    n_gcpcs <- inspected$n_gcpcs
    
    # Compute covariance matrices
    RaRa <- crossprod(Ra) / (nrow(Ra) - 1)
    RbRb <- crossprod(Rb) / (nrow(Rb) - 1)
    
    # --- Simple contrastive version ---
    if (method == 'v1') {
      JRaRaJ <- t(J) %*% RaRa %*% J
      JRbRbJ <- t(J) %*% RbRb %*% J
      sigma <- JRaRaJ - alpha * JRbRbJ
      eig <- eigen(sigma, symmetric = TRUE)
      w <- eig$values
      v <- eig$vectors
      eig_idx <- order(w, decreasing = TRUE)
      x <- J %*% v[, eig_idx]
      s_total <- w[eig_idx]
      obj_info <- 'Ra - alpha * Rb'
      
      # --- Generalized contrastive PCA versions ---      
    } else {
      denom_well_conditioned <- FALSE
      ortho_column_order <- c()
      count_dim <- 0
      x <- NULL
      x_orth <- NULL
      
      # Iterate to extract components one by one
      for (idx in 1:n_gcpcs) {
        JRaRaJ <- t(J) %*% RaRa %*% J
        JRbRbJ <- t(J) %*% RbRb %*% J
        
        # Define numerator and denominator for generalized Rayleigh quotient
        if (method %in% c('v2', 'v2.1')) {
          numerator <- JRaRaJ; denominator <- JRbRbJ
          obj_info <- 'Ra / Rb'
        } else if (method %in% c('v3', 'v3.1')) {
          numerator <- JRaRaJ - JRbRbJ; denominator <- JRbRbJ
          obj_info <- '(Ra-Rb) / Rb'
        } else if (method %in% c('v4', 'v4.1')) {
          numerator <- JRaRaJ - JRbRbJ; denominator <- JRaRaJ + JRbRbJ
          obj_info <- '(Ra-Rb) / (Ra+Rb)'
        } else stop("Unrecognized method.")
        
        if (!denom_well_conditioned) {
          if (kappa(denominator) > cond_number) {
            warning("Denominator ill-conditioned. Regularizing.")
            w <- eigen(denominator, symmetric = TRUE)$values
            alpha <- max(w[1] / cond_number - tail(w, 1), 0)
            denominator <- denominator + diag(nrow(denominator)) * alpha
          }
          denom_well_conditioned <- TRUE
        }
        
        # Whitening the denominator
        eig_den <- eigen(denominator, symmetric = TRUE)
        d <- eig_den$values
        d[d < 0] <- 0
        e <- eig_den$vectors
        M <- e %*% diag(sqrt(d)) %*% t(e)
        Minv <- tryCatch(solve(M), error = function(e) solve(M + diag(nrow(M)) * 1e-6))
        sigma <- t(Minv) %*% numerator %*% Minv
        
        # Solve eigenproblem in transformed space
        eig <- eigen(sigma, symmetric = TRUE)
        w <- eig$values
        v <- eig$vectors
        eig_idx <- order(w, decreasing = TRUE)
        v <- v[, eig_idx]
        x_temp <- J %*% Minv %*% v
        x_temp <- x_temp / norm(x_temp, type = "2")
        
        # Orthogonalize and select directions
        if (idx == 1) {
          x <- x_temp
          x_orth <- matrix(x_temp[, 1], ncol = 1)
          ortho_column_order <- c(ortho_column_order, count_dim)
          count_dim <- count_dim + 1
        } else {
          if (idx %% 2 == 1) {
            x_add <- matrix(x_temp[, ncol(x_temp)], ncol = 1)
            ortho_column_order <- c(ortho_column_order, ncol(x_temp) + count_dim - 1)
          } else {
            x_add <- matrix(x_temp[, 1], ncol = 1)
            ortho_column_order <- c(ortho_column_order, count_dim)
            count_dim <- count_dim + 1
          }
          x_orth <- cbind(x_orth, x_add)
        }
        
        if ((ncol(J) - ncol(x_orth)) < 1) break  # avoid shrinking to zero cols
        j <- svd(J - x_orth %*% t(x_orth) %*% J)$u
        J <- j[, 1:(ncol(J) - 1), drop = FALSE]
      }
      
      if (method %in% c('v2.1', 'v3.1', 'v4.1')) {
        new_column_order <- order(ortho_column_order)
        x <- x_orth[, new_column_order, drop = FALSE]
      }
      
      # Compute scores and objective values
      RaX <- Ra %*% x
      RbX <- Rb %*% x
      XRaRaX <- crossprod(RaX)
      XRbRbX <- crossprod(RbX)
      
      if (method %in% c('v2', 'v2.1')) {
        numerator_orig <- XRaRaX; denominator_orig <- XRbRbX
      } else if (method %in% c('v3', 'v3.1')) {
        numerator_orig <- XRaRaX - XRbRbX; denominator_orig <- XRbRbX
      } else if (method %in% c('v4', 'v4.1')) {
        numerator_orig <- XRaRaX - XRbRbX; denominator_orig <- XRaRaX + XRbRbX
      }
      
      s_total <- diag(numerator_orig) / diag(denominator_orig)
    }
    
    # Return loadings, scores, and diagnostics
    list(loadings = x,
         Ra_scores = Ra %*% x / norm(Ra %*% x, type = "2"),
         Ra_values = norm(Ra %*% x, type = "2"),
         Rb_scores = Rb %*% x / norm(Rb %*% x, type = "2"),
         Rb_values = norm(Rb %*% x, type = "2"),
         objective_function = obj_info,
         objective_values = s_total,
         null_gcpca_values = if (Nshuffle > 0) null_distribution(Ra, Rb) else NULL)
  }
  
  null_distribution <- function(Ra, Rb) {
    null_vals <- NULL
    for (i in 1:Nshuffle) {
      Ra_shuf <- Ra[sample(nrow(Ra)), ]
      Rb_shuf <- Rb[sample(nrow(Rb)), ]
      fit_result <- fit(Ra_shuf, Rb_shuf)
      null_vals <- rbind(null_vals, fit_result$objective_values)
    }
    null_vals
  }
  
  fit(Ra, Rb)
}