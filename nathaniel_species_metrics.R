library(dplyr)

taxon_group <- "eucalypts_2016Phylo"
taxon_dir <- paste0("Z:/!my_working/data/", taxon_group, "/")
data_file <- paste0(taxon_dir, taxon_group, "_almalgamated_data.csv")

data <- read.csv(data_file, header=T,sep=",")



species_list <- data[!grepl("\\d",data$node), ]


# View(species_list)

map_equation_conv_flow <- arrange(species_list, desc(map_equation_conv_flow))
modularity_conv_participation_coefficient <- arrange(species_list, desc(modularity_conv_participation_coefficient))
modularity_conv_within_module_degree <- arrange(species_list, desc(modularity_conv_within_module_degree))

network_conv_betweenness <- arrange(species_list, desc(network_conv_betweenness))
network_rw_betweenness <- arrange(species_list, desc(network_rw_betweenness))

data_subset <- network_rw_betweenness
head(data_subset$node, 3)
tail(data_subset$node, 3)
