### Get traits with the strongest effect sizes.

pos3_list <- lapply(cohen_list, function(df) {
  df <- df[!(rownames(df) %in% total_traits),] # Remove 'severity_score' and 'total' traits
  rownames(df)[order(as.numeric(df[,2]), decreasing = TRUE)[1:5]] # Get top 3 positive values
})

neg3_list <- lapply(cohen_list, function(df) {
  df <- df[!(rownames(df) %in% total_traits),]
  rownames(df)[order(as.numeric(df[,2]), decreasing = FALSE)[1:5]]
})

# Assign names to the lists
names(pos3_list) <- disease_name
names(neg3_list) <- disease_name

# Get unique values across all lists
pos3 <- unique(unlist(pos3_list))
neg3 <- unique(unlist(neg3_list))


# Create a new column 'family' for disease family in the dataframe
cohen_long$family <- NA
for (family_name in names(families)) {
  cohen_long$family[cohen_long$trait %in% families[[family_name]]] <- family_name
}

# Create a new column 'type' for trait family in the dataframe
cohen_long$type <- NA
for (disease_type in names(disease_types)) {
  cohen_long$type[cohen_long$disease %in% disease_types[[disease_type]]] <- disease_type
}

# Reorder the factor levels to adjust legend order
cohen_long$disease <- factor(cohen_long$disease, levels = disease_order)


# Keep only the strongest effect sizes
subset_pos_cohen <- do.call(rbind, lapply(names(pos3_list), function(disease_name) {
  subset(cohen_long, disease == disease_name & trait %in% pos3_list[[disease_name]])
}))

subset_neg_cohen <- do.call(rbind, lapply(names(neg3_list), function(disease_name) {
  subset(cohen_long, disease == disease_name & trait %in% neg3_list[[disease_name]])
}))

subset_cohen <- rbind(subset_pos_cohen, subset_neg_cohen)