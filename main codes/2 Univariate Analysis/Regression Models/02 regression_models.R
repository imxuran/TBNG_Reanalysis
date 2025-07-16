### Regresssion models

# Create new dataframes to store the results of three models
col_names <- c("trait", "estimate", "SE", "stas", "p_value", "CI_lower", "CI_upper")
age_df <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(age_df) <- col_names
sex_df <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(sex_df) <- col_names
both_df <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(both_df) <- col_names



for (t in common_traits) {
  
  # Linear model for age
  lm_age <- lm(as.formula(paste("age ~", t)), data = df_std)
  age_df[nrow(age_df) + 1, ] <- tidy(lm_age, conf.int = TRUE, exponentiate = FALSE) %>%
    filter(term != "(Intercept)")
  
  
  # Logistic regression model for sex
  lm_sex <- glm(as.formula(paste("sex ~", t)), family = binomial, data = df_std)
  sex_df[nrow(sex_df) + 1, ] <- tidy(lm_sex, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)")
  
  
  # Linear models with interaction terms
  lm_both <- lm(as.formula(paste("age ~", t, "* sex")), data = df_std)
  both_df[nrow(both_df) + 1, ] <- tidy(lm_both, conf.int = TRUE, exponentiate = FALSE) [4,]
  
}


# Reorder the daraframes based on the absolute value of the estimates
age_df <- age_df[order(abs(age_df$estimate), decreasing = TRUE), ]
sex_df <- sex_df[order(abs(sex_df$estimate), decreasing = TRUE), ]
both_df <- both_df[order(abs(both_df$estimate), decreasing = TRUE), ]
both_df$trait<-sub(":sex$", "", both_df$trait)