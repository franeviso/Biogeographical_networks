##########################################################################################################
# 
#
#Nunzio.Knerr@csiro.au
#updated Date:30/10/2015
#
# Edited by Nathaniel Bloomfield
##########################################################################################################
library(tools)


# inputs
#taxon_group <- "acaciaPhylo" 
#taxon_folder <- paste0("Z:/!my_working/data/", taxon_group)

#analysis_type <- "map" #map or modularity
#analysis <- "phylo" #conventional, phylo or phylo_summed

args <- commandArgs(trailingOnly = TRUE) 
taxon_group <- args[1] 
analysis_type <- args[2] #map or modularity
analysis <- args[3]  #conventional, phylo or phylo_summed
iterations <- as.numeric(args[4]) 

taxon_folder <- paste0("~/phylogenetic_network/", taxon_group)
# map equation...

map_equation_install_folder <- paste0("~/phylogenetic_network/map_equation/Infomap/")
number_of_iterations <- paste0("-N ", iterations, " ")
other_parameters <- paste0("--undirected --two-level ")
output_formats <- paste0("--tree --bftree --map --clu ")
map_arguments <- paste0(number_of_iterations, other_parameters, output_formats)

# modularity equation...
modularity_install_folder <- paste0("~/phylogenetic_network/netcarto_cl/")
T_ini <- paste0(5, " ") #2/size_of_network default
iteration_factor <- paste0(1, " ") #range from 0.1, 1
cooling_factor <- paste0(0.95, " ") #min 0.9 and max 0.995
randomizations <- iterations
mod_arguments <- paste0(T_ini, iteration_factor, cooling_factor, randomizations)


# do not modify #########################################

# Functions

InvokeMap <- function(analysis,taxon_group, taxon_folder, map_arguments, map_equation_install_folder) {
  
  if(analysis == "phylo_summed"){
    suffix <- "_all_branch_pajek_summed_branch_length.net"
    output_name <- "output_map_phylo_summed"
    input_name <- "for_network_analyses"
  } else if (analysis == "phylo") {
    suffix <- "_all_branch_pajek_single_branch_length.net"
    output_name <- "output_map_phylo"
    input_name <- "for_network_analyses"
  } else if (analysis == "conventional") {
    suffix <- "_tips_only_pajek_no_weights.net"
    output_name <- "output_map"
    input_name <- "for_network_analyses_tips_only"
  }
  
  input_file <- paste0(taxon_folder,"/",input_name, "/", taxon_group, suffix)
  output_folder <- paste0(taxon_folder, "/", output_name, '/')
  dir.create(output_folder, showWarnings = FALSE)
  
  seed <- paste0("--seed ", sample(1:10000, 1), " ")  # random number between 1 and 10000
  
  ptm <- proc.time()#start timer
  
  cmd <- paste(map_equation_install_folder, "Infomap ", input_file, " ", output_folder, " ", seed, map_arguments, sep="")
  system(cmd, intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE) 
  
  sink(paste(taxon_folder, taxon_group, output_name, taxon_group, "_network_analysis_run_stats.txt", sep=""))
  cat("Summed branch length analyses (map equation): ")
  cat("\n")
  cat(cmd)
  cat("\n")
  cat("analysis took: ")
  cat(proc.time() - ptm) #stop timer
  cat("\n")
  cat(paste0("results should be in: ", output_folder))
  cat("\n")
  cat("\n")
  cat("##########################################")
  cat("\n")
  sink()
  
}


InvokeModularity <- function(analysis,taxon_group, taxon_folder, mod_arguments, modularity_install_folder) {
  
  
  if(analysis == "phylo_summed" || analysis == "phylo" ){
    suffix <- "_all_branch_links_list_bipartite_network.txt"
    output_name <- "output_mod_phylo"
    input_name <- "for_network_analyses"
  } else if (analysis == "conventional") {
    suffix <- "_tips_only_links_list_bipartite_network.txt"
    output_name <- "output_mod"
    input_name <- "for_network_analyses_tips_only"
  }
  input_file <- paste0(taxon_folder,"/",input_name, "/", taxon_group, suffix)
  
  output_folder <- paste0(taxon_folder, "/", output_name, '/') # the directory that the results will go to
  dir.create(output_folder, showWarnings = FALSE)
  
  seed <- paste0(sample(1:10000, 1), " ")  # random number between 1 and 10000

  
  ptm <- proc.time()#start timer
  
  setwd(output_folder)
  file_prefix <- basename(file_path_sans_ext(input_file))#file name only with no extension
  
  cmd <- paste(modularity_install_folder, "netcarto_cl ", input_file, " ", seed, mod_arguments, sep="")
  system(cmd, intern=FALSE, ignore.stdout=FALSE, ignore.stderr=FALSE, wait=TRUE) 
  
  #getwd()
  file_list <- dir(getwd())
  file.rename(file_list, paste0(file_prefix, " ", file_list))
  
  sink(paste(taxon_folder, taxon_group, output_name, taxon_group, "network_analysis_run_stats.txt", sep=""), append=TRUE)
  cat("\n")
  cat("conventional occurrence network analyses (modularity): ")
  cat("\n")
  cat(cmd)
  cat("\n")
  cat("analyses took: ")
  cat(proc.time() - ptm)#stop timer
  cat("\n")
  cat(paste0("results should be in: ", output_folder))
  cat("\n")
  cat("\n")
  cat("##############################################")
  cat("\n")
  sink()
}

# Outputs

if (analysis_type == 'map') {
  InvokeMap(analysis, taxon_group, taxon_folder, map_arguments, map_equation_install_folder)
} else if (analysis_type == 'modularity') {
  InvokeModularity(analysis, taxon_group, taxon_folder, mod_arguments, modularity_install_folder)
}





