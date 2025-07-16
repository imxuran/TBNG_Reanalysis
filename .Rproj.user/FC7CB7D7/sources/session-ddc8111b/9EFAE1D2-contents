### Calculate the effect size, obtain one list including information on cohen'd, SE and CI for each dataset.


# Function to calculate the baised std
biased_sd <- function(x) {
  sqrt(sum((x - mean(x))^2) / length(x))
}

# Function to calculate Cohen's d and related statistics
calculate_cohens_d <- function(df) {
  
  # Select numeric columns excluding age
  numeric_cols <- sapply(df, is.numeric) & !(tolower(names(df)) %in% c("age","weight","height","pdac_binary","severity_score","sample id"))
  
  # Get two groups, only save the traits
  case <- df[df$group=="case", numeric_cols]
  #intersect(colnames(df), setdiff(trait_order,total_traits))]
  control <- df[df$group=="control", numeric_cols]
  #intersect(colnames(df), setdiff(trait_order,total_traits))]
  
  # Calculate means, biased standard deviations, and sample sizes
  mean_case <- apply(case, 2, mean)
  mean_control <- apply(control, 2, mean)
  sd_case <- apply(case, 2, biased_sd)
  sd_control <- apply(control, 2, biased_sd)
  len_case <- nrow(case)
  len_control <- nrow(control)
  
  # Calculate pooledSD, Cohen's d, standard error, and confidence intervals
  poolSd <- sqrt( ((len_case-1) * sd_case^2 + (len_control-1) * sd_control^2)  / (len_case+len_control-2) )
  cohenD <- (mean_case - mean_control) / poolSd
  SE_d <- sqrt( (len_case + len_control)/(len_case * len_control) + cohenD^2/(2*(len_case+len_control)) )
  CI_lower <- cohenD - 1.96 * SE_d
  CI_upper <- cohenD + 1.96 * SE_d
  
  # Return a data frame for each disease with traits as rownames and stats as columns
  data.frame(
    trait = colnames(case),
    cohen_d = cohenD,
    SE = SE_d,
    CI_lower = CI_lower,
    CI_upper = CI_upper
  )
}


# List to store the data frames for each disease
cohen_list <- lapply(data_list, calculate_cohens_d)

# Reshape the datasets into a dataframe compatible with ggforestplot.
cohen_long <- bind_rows(cohen_list, .id = "disease")