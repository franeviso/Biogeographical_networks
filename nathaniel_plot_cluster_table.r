library(dplyr)

taxon_group <- "eucalypts_2016Phylo"
# NAME <- 'map_equation_conv_moduleID'
# Sorter <- 'map_equation_conv_flow'
# TableTitle <- "Map Equation Bioregions"

taxon_dir <- paste0("Z:/!my_working/data/", taxon_group, "/")

for_rast <- paste0(taxon_dir, taxon_group, "_almalgamated_data_geo_only.csv")
for_rast_sum <- paste0(taxon_dir, taxon_group, "_almalgamated_data.csv")
summary_data <- read.csv(for_rast_sum, header=T,sep=",")


output_folder <- paste0(taxon_dir, "tables")
dir.create(output_folder, showWarnings = FALSE)


species_list <- summary_data[!grepl("\\d",summary_data$node), ]

locations_list <- read.csv(for_rast, header=T,sep=",")
# locations_list <- summary_data[grepl("p_\\d",summary_data$node), ]
# View(locations_list)

for (i in 1:3) {
  if (i == 1) {
    NAME <- 'map_equation_conv_moduleID'
    Sorter <- 'map_equation_conv_flow'
    TableTitle <- "Map Equation Bioregions. Indicative species are those with the highest
     flow. The codelength for this partition was X."
  } else if (i == 2) {
    NAME <- 'map_equation_rw_moduleID'
    Sorter <- 'map_equation_rw_flow'
    TableTitle <- "Map Equation RW Bioregions. Indicative species are those with the highest
     RW flow. The codelength for this partition was X."
    
  } else if (i == 3) {
    NAME <- 'modularity_conv_moduleID'
    Sorter <- 'modularity_conv_within_module_degree'
    TableTitle <- "Modularity Bioregions. Indicative species are those with the highest
     within module degree. The randomisation test for this network found that the average
     modularity was X, with a standard deviation of X. The modularity value for this 
     partition was X, which had a p value of X"
    
  }
  regions <- arrange_(distinct_(select_(species_list, NAME)), NAME)
  
  output_summary <- paste0(taxon_group, '_', NAME, "_summary.txt")
  sink(paste0(output_folder,'/', output_summary))
  
  cat("\\begin {table}[h]")
  cat(paste0("\\caption {", TableTitle,"} \n")) 
  cat("\\begin{center} \n")
  cat("\\footnotesize \n")
  cat("\\begin{tabular}{ |lllll| } \n")
  cat("\\hline \n")
  
  # cat(" Region No. & No. of Species & No. of Locations & Relative Diversity & Sum Inverse Range & Avg. Inverse Range & Ind. Species \\\\ \n")
  
  cat(" Region No. & No. Species & No. Locations & Avg. Inv. Range & Ind. Species \\\\ \n")
  cat("\\hline \n")
  
  for (region in 1:length(regions[[NAME]])) {
    in_region <- as.data.frame(species_list[species_list[, NAME] == toString(regions[[NAME]][region]), ])
    in_region_locations <- as.data.frame(locations_list[locations_list[, NAME] == toString(regions[[NAME]][region]), ])
    #         if (region == 1) {
    #           View(in_region_locations)
    #         }
    indicative_species <- head(in_region[order(desc(in_region[[Sorter]])), ], 2) 
    
    region_WE <- sum(1/in_region$modularity_conv_num_links)
    region_CWE <- region_WE/nrow(in_region)
    #         region_WE <- sum(in_region$network_conv_degree)
    #         region_CWE <- 1/region_WE/nrow(in_region)
    sig_figs <- 3
  #   cat(paste0( region, " & ", nrow(in_region), " & ", nrow(in_region_locations), " & ", signif(nrow(in_region)/nrow(in_region_locations), digits = sig_figs),
  #               " & ",  signif(region_WE, digits = sig_figs)," & ",signif(region_CWE, digits = sig_figs)," & ",noquote(toString(indicative_species$node)), "\\\\ \n"))
    table_row <- paste0( region, " & ", nrow(in_region), " & ", nrow(in_region_locations), " & ", 
                signif(region_CWE, digits = sig_figs)," & ",noquote(toString(indicative_species$node)), "\\\\ \n")
    table_row <- gsub("_", " ", table_row, fixed = TRUE)
    cat(table_row)
    cat("\\hline \n")
  }
  
  cat("\\end{tabular} \n")
  cat("\\end{center} \n")
  cat("\\end {table}")
  sink()
}
# View(new_all_cells)