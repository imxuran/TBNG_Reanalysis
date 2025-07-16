### DMFA implementation and results extraction

# Merge data as input
df_all_case <- bind_rows(
  lapply(data_list, function(df) df[df$group=="case", analy_traits]), .id = "disease")

df_all <- bind_rows(
  lapply(data_list, function(df) df[, analy_traits]), .id = "disease")

# conduct DMFA on cases
dmfa_case <- DMFA(df_all_case, 
                  num.fact = 1,
                  graph = FALSE)

# conduct DMFA on controls
dmfa_all <- DMFA(df_all, 
                 num.fact = 1,
                 graph = FALSE)


# Extract Variance
eigenvalues <- dmfa_all$eig[,1]
prop_var <- eigenvalues / sum(eigenvalues) * 100


## Score data
# Convert DMFA coordinates to a dataframe
dmfa_coords <- data.frame(
  x = dmfa_all$ind$coord[,1],
  y = dmfa_all$ind$coord[,2],
  group = df_all$disease
)

# Change the disease order
dmfa_coords$group <- factor(dmfa_coords$group, levels = disease_order)


## Loading data
# Combine data from all groups
dmfa_loadings <- data.frame()
for(i in 1:length(dmfa_all$var.partiel)) {
  loadings <- as.matrix(dmfa_all$var.partiel[[i]])
  group_data <- data.frame(
    Variable = rownames(loadings),
    Dim1 = loadings[, 1],
    Dim2 = loadings[, 2],
    Group = names(dmfa_all$var.partiel)[[i]]
  )
  dmfa_loadings <- rbind(dmfa_loadings, group_data)
}

# Create unit circle coordinates
theta <- seq(0, 2*pi, length.out = 100)
circle_data <- data.frame(
  x = cos(theta),
  y = sin(theta)
)

dmfa_loadings$Group <- factor(dmfa_loadings$Group, levels = unlist(disease_types))


