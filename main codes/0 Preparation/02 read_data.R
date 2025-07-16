### Read data from the folders, the raw data should be put in the 'data/data_part_II' filecase,
### while the derived data should be put in the 'data/TPNG_data_part_I' filecase. 


# Path of the derived and original data set
files <- list.files(path = "data/TBNG_data_part_I", pattern = "\\.rds$", full.names = TRUE)
files_direct <- list.files(path = "data/TBNG_data_part_II", pattern = "\\.rds$", full.names = TRUE)


# Read data
data_list <- lapply(files, readRDS)
data_list_direct <- lapply(files_direct, readRDS)


# Extract disease names for derived and original set
disease_name <- sapply(basename(files), function(file) {
  tail(strsplit(tools::file_path_sans_ext(file), "_")[[1]], 1)
})

disease_name_2 <- sapply(basename(files_direct), function(file) {
  strsplit(tools::file_path_sans_ext(file), "_")[[1]][[2]]
})


# Assign name to every data frame stored in the list
names(data_list) <- disease_name
names(data_list_direct) <- unique(disease_name_2)
