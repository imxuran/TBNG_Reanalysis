### Define the variables used multiple times in analysis



## 1 TYPES

# Define the diseases from different types
disease_types <- list(
  C = c("MM", "CRC", "PC", "PDAC", "BC"),
  AI = c("AIH","CD","UC","RA"),
  INF = c("LMS","COVID19"),
  M = c("T2D", "MASH")
)

# Define the derived traits from different Families
families <- list(
  Complexity = c("CA1", "CA2", "CA3", "CA4", "TC", "TM", "THy", "MM"),
  Bisection = c("A2B", "CB", "TB", "A2FS0B"),
  Galactosylation = c("A2FG", "A2FS0G", "A2G", "CG","A4G"),
  Fucosylation = c("A1F", "A1Fa", "A2F", "A2Fa", "A3F", "A3Fa", "A4F", "A4Fa", "CF", "CFa", "TF"),
  Sialylation = c("A2S", "A3S", "A4S", "CS", "TA1S", "TA2S", "TA3S", "TA4S"),
  Sialylation_L = c("A1L", "A2L", "A3L", "A4L"),
  Sialylation_E =c("A1E", "A2E", "A3E", "A4E")
)



## 2 COLORS

# Define the colors representing diseases
disease_colors <- c(AIH = "#CBC9DA", RA = "#9D98BE", CD = "#6C64A2", UC = "#4C3E93", COVID19 = "#EA9C9D", LMS = "#DB6968", T2D = "#99CBEB", MASH = "#4D97CD", CRC = "#7D9F84", PC = "#93CC82", PDAC = "#74C476", MM = "#459943", BC = "#3B754A")

# Define the colors representing disease types
type_colors <- c(AI = "#9D98BE", INF = "#EA9C9D", M = "#99CBEB", C = "#93CC82")

# Define the colors of different families
family_colors <- c(Complexity = "#6BCCEF90", Bisection = "#00D187", Galactosylation = "#009462", Fucosylation = "#F6AF03", Sialylation = "#9C9C9C", Sialylation_L = "#951255", Sialylation_E = "#DB6968")



## 3 ORDERS

# Specify the new order of traits and diseases, make traits/diseases from the same family/type together
trait_order <- unlist(families)
disease_order <- unlist(disease_types)



## 4 TRAITS

# All present traits in 9 original datasets, filter out basic information, only keep traits whose name starting with 'H'
all_direct_traits <- unique(grep("^H", unlist(lapply(data_list_direct, colnames)), value = TRUE))

# Common derived traits in all datasets
common_traits <- Reduce(intersect, lapply(data_list, names))
common_traits <- common_traits[common_traits != "group"]

# Exclude the "total" derived traits
total_traits <- c("TA1S", "TA2S", "TA3S", "TA4S", "TA1", "TA2", "TA3", "TA4", "TF", "TB")

# Traits used in analysis (common traits excluding 'total' traits)
analy_traits <- setdiff(common_traits,total_traits)