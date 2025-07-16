### Calculate variable similarity matrices per dataset (case and control)


## For case samples

# Only keep case group
case_list <- lapply(data_list, function(df) {
  df<-df[df$group=='case', analy_traits]
})

# list to keep similarity matrix of each dataset
similarity_case <- list()

# Pearson used in the calculation
for (df in case_list) {
  
  s <- compute_similarity(df)
  
  similarity_case <- append(similarity_case, list(s))
}

names(similarity_case) <- disease_name


## For control samples

# Only keep control group
control_list <- lapply(data_list, function(df) {
  df<-df[df$group=='control', analy_traits]
})

# list to keep similarity matrix of each dataset
similarity_control <- list()

# Pearson used in the calculation
for (df in control_list) {
  
  s <- compute_similarity(df)
  
  similarity_control <- append(similarity_control, list(s))
}

names(similarity_control) <- disease_name