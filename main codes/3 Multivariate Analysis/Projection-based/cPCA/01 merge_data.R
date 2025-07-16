### Prepare data


# Merge all data together, grouped by disease status
df_all_case <- bind_rows(
  lapply(data_list, function(df) df[df$group=="case", analy_traits]), .id = "disease")

df_all_control <- bind_rows(
  lapply(data_list, function(df) df[df$group=="control", analy_traits]), .id = "disease")

df_all_control$disease <- 'control'

df_all <- rbind(df_all_case,df_all_control)

group <- df_all$disease

data_numeric <- df_all[,-1]
