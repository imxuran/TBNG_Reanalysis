### Use Upset plot to compare and assess the measured glycan compositions across datasets involved.


## Create a binary matrix to indicate trait presence, 
## check which traits from all_traits are present in the current disease dataframe.
trait_presence_matrix <- sapply(data_list_direct, function(x) {
  as.integer(all_direct_traits %in% colnames(x))
})

trait_presence_df <- as.data.frame(trait_presence_matrix)

rownames(trait_presence_df) <- all_direct_traits


## Define metadata for further plotting
sets <- unique(ifelse(disease_order %in% c('CD', 'UC'), 'IBD', disease_order))
types <- rep(names(disease_types), sapply(disease_types, length))
types <- types[-match("AI", types)]
metadata <- as.data.frame(cbind(sets, types))
names(metadata) <- c("sets", "disease")


## Plot the upset plot

png("Results/Compositions_Upset.png", width = 5500, height = 4000, res = 500)

upset(trait_presence_df, nintersects = 10, 
      point.size = 2.7, line.size = 1, 
      mainbar.y.label = "Intersections",
      sets.x.label = "Glycan Compositions", 
      text.scale = c(1.5, 1.5, 1.3, 1.5, 1.5, 1.5), 
      matrix.color = "#0074B3", sets.bar.color = "#F47720", main.bar.color = "#F47720",
      keep.order = TRUE, order.by = "freq", sets = sets,
      #sets = c("AIH", "IBD", "RA", "LMS", "COVID19", "T2D", "MASH", "BC", "CRC", "MM", "PC"),
      set.metadata = list(data = metadata, plots = list(
        list(type = "text", column = "disease", assign = 3, colors = type_colors, title = NULL),
        list(type = "heat", column = "disease", assign = 7, colors = type_colors, alpha = 0.5, title = NULL),
        list(type = "matrix_rows", column = "disease", colors = type_colors, alpha = 0.5))),
      queries = list(list(query = intersects, params = sets, color="#0074B3",  active = T))
)

dev.off()