### Pool all control dataset together


data_list <- lapply(data_list, function(df) {
  # Rename all age columns as "age", some are 'Age' / 'AGE'
  df <- rename_with(df, tolower, .cols = matches("^age$", ignore.case = TRUE))
  names(df)[str_detect(names(df), regex("sex", ignore_case = TRUE))] <- "sex"
  # Convert a factor to numeric if 'age' exsits
  if ("age" %in% colnames(df)) {
    df$age <- as.numeric(as.character(df$age))
  }
  return(df)
})


# Combine all dataframes
df <- bind_rows(data_list)


# Only keep control group, column: age, sex, common traits
df <- df[df$group=='control', c("age","sex", common_traits)] %>%
  # Standardize sex column
  mutate(sex = case_when(
    sex %in% c('female', 'F') ~ 1,
    sex %in% c('male', 'M') ~ 0,
    TRUE ~ NA_real_
  ))


# Remove rows without age or sex
df <- na.omit(df)

# Standardize all traits
df_std <- cbind(df[,1:2],scale(df[,-(1:2)]))