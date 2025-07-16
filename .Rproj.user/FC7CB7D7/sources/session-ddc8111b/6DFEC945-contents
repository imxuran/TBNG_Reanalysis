### WGCNA implementation and UpSet plot for assessing the most disease-associated modules for each study


# implement WGCNA and get the most disease-associated traits
module_genes <- lapply(names(data_list), function(n) wgcna(df = data_list[[n]], n = n))

names(module_genes) <- names(data_list)

# Integrate the traits' presence
module_presence_matrix <- sapply(module_genes, function(x) {
  as.integer(common_traits %in% x)
})

module_presence_df <- as.data.frame(module_presence_matrix)

rownames(module_presence_df) <- common_traits

# Metadata for upset plot annotation
types <- rep(names(disease_types), sapply(disease_types, length))
metadata <- as.data.frame(cbind(disease_order, types))
names(metadata) <- c("sets", "disease")

png("Results/WGCNA/individual_Upset.png", width = 5500, height = 4000, res = 500)

upset(module_presence_df, 
      #nintersects = 10, 
      point.size = 2.7, line.size = 1, 
      mainbar.y.label = "Intersections",
      sets.x.label = "Glycan Compositions", 
      text.scale = c(1.5, 1.5, 1.3, 1.5, 1.5, 1.5), 
      matrix.color = "#0074B3", sets.bar.color = "#F47720", main.bar.color = "#F47720",
      keep.order = TRUE, order.by = "degree", show.numbers = FALSE,
      #sets = sets,
      sets = disease_order, set_size.scale_max = 13,
      set.metadata = list(data = metadata, plots = list(
        list(type = "text", column = "disease", assign = 3, colors = type_colors, title = NULL),
        list(type = "heat", column = "disease", assign = 7, colors = type_colors, alpha = 0.5, title = NULL),
        list(type = "matrix_rows", column = "disease", colors = type_colors, alpha = 0.5))),
      queries = list(
        list(query = intersects, params = c('MASH','UC','PDAC','MM'), color="#951255",  active = T), # A2L
        list(query = intersects, params = c('T2D','COVID19','RA','CD','AIH','BC','PC','MM'), color="#9C9C9C",  active = T), #A2E+CS
        list(query = intersects, params = c('T2D','RA','CD','AIH','BC','CRC'), color="#00D187",  active = T), #Bisection
        list(query = intersects, params = c('T2D','RA','CD','AIH','BC','CRC','MM'), color="#F6AF03",  active = T), # Bi-Fuco
        list(query = intersects, params = c('T2D','COVID19','CD','AIH','BC','PC','MM'), color="#009462",  active = T), #Bi-Gala
        list(query = intersects, params = c('MASH','LMS','UC','PDAC'), color="#DB6968",  active = T), #A3E
        list(query = intersects, params = c('MASH','UC','PDAC'), color="#9C9C9C",  active = T), #Other Sia
        list(query = intersects, params = c('LMS','CD','CRC'), color="#6BCCEF",  active = T), #Complex T
        list(query = intersects, params = c('LMS','CD','AIH'), color="#6BCCEF",  active = T), #MM
        list(query = intersects, params = c('UC'), color="#F6AF03",  active = T) #Other Fuco
      )
)


dev.off()