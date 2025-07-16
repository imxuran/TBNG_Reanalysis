### Calculate variances of coordinates for variables

dims = c(1, 2)
n_groups <- length(dmfa_case$var.partiel)
n_vars <- nrow(dmfa_case$var.partiel[[1]])

# Store coordinates for each variable across groups
coords_by_var <- array(NA, dim = c(n_vars, 2, n_groups))

# Extract coordinates for each group
for(i in 1:n_groups) {
  coords <- as.matrix(dmfa_case$var.partiel[[i]])
  coords_by_var[,,i] <- coords[, dims]
}

# Calculate variance in position for each variable
position_variance <- data.frame(
  Variable = rownames(dmfa_case$var.partiel[[1]]),
  Variance_Dim1 = apply(coords_by_var[,1,], 1, var),
  Variance_Dim2 = apply(coords_by_var[,2,], 1, var),
  Total_Variance = apply(coords_by_var[,1,], 1, var) + 
    apply(coords_by_var[,2,], 1, var)
)

# Sort by total variance
position_variance <- position_variance[order(position_variance$Total_Variance, decreasing = TRUE),]

print(position_variance)