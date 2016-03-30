library(tools)
library(blockmodeling)

taxon_group <- "acaciaPhylo"

taxon_dir <- paste0("Z:/!my_working/data/", taxon_group, "/")

# locations_only <- T
do_networks <- T
do_biodiverse <- T

conventional <- T
phylo <- T
mod_phylo <- T
range_weighted <- T
range_weighted_phylo <- T

mod <- T
map <- T
########################################################
# load the map equation output
########################################################

Map_analysis <- function(wd, taxon_group, suffix, name) {
  setwd(wd)
  
  results_clu <- paste0(taxon_group, name,".clu")
  results_tree <- paste0(taxon_group, name, ".tree")
  
  map_clu <- loadvector(results_clu)
  colnames(map_clu) <- c("ID",	paste0("map_equation_", suffix, "_moduleID"),	paste0("map_equation_", suffix, "_flow"))

  # load the map equation names and ID
  map_nameID <- read.table(results_tree, sep = "" , header = F, skip = 1)
  map_nameID <- map_nameID[,-(1:2)] #
  colnames(map_nameID) <- c("node",	"ID")
  total_map <- Reduce(function(...) merge(..., all=T), list(map_nameID, map_clu))
  return(total_map)
}


total_map <- data.frame(node = '')
total_map_phylo <- data.frame(node = '')
total_map_rw <- data.frame(node = '')
total_map_prw <- data.frame(node = '')
if (map == T) {
  if (conventional == T) {
    total_map <- Map_analysis(paste0(taxon_dir, "output_map/"), taxon_group, 'conv', "_tips_only_pajek_no_weights")
  } 
  if (phylo == T) {
    total_map_phylo <- Map_analysis(paste0(taxon_dir, "output_map_phylo/"), taxon_group, 'phylo', "_all_branch_pajek_single_branch_length")
  }
  if (range_weighted == T) {
    total_map_rw <- Map_analysis(paste0(taxon_dir, "output_map_range_weighted/"), taxon_group, 'rw', "_tips_only_pajek_range_weighted")
  }
  if (range_weighted_phylo == T) {
    total_map_prw <- Map_analysis(paste0(taxon_dir, "output_map_phylo_range_weighted/"), taxon_group, 'prw', "_all_branch_pajek_single_branch_length_range_weighted")
  }
}

##################################################
# load the modularity output
########################################################

Mod_analysis <- function(wd, taxon_group, suffix, name) {
  setwd(wd)
  results_modules <- paste0(taxon_group, name, "modules.clu")
  results_roles <- paste0(taxon_group, name, "roles.clu")
  results_prop <- paste0(taxon_group, name, "node_prop.dat")
  
  
  # grabbing modularity moduleID and ID
  mod_clu <- read.table(results_modules, sep = "" , header = F, skip = 1)
  colnames(mod_clu) <- c(paste0("modularity_", suffix,"_moduleID"))
  mod_clu$ID<-seq.int(nrow(mod_clu))
  
  # grabbing the role of the node in the network
  mod_clu_role <- read.table(results_roles, sep = "" , header = F, skip = 1)
  colnames(mod_clu_role) <- c(paste0("modularity_", suffix, "_role"))
  mod_clu_role$ID<-seq.int(nrow(mod_clu_role))
  
  
  # grabbing the number of links, within module degree and participation coefficient
  mod_clu_prop <- read.table(results_prop, sep = "" , header = F, skip = 0)
  colnames(mod_clu_prop) <- c("node", paste0("modularity_", suffix, "_num_links"), paste0("modularity_", suffix, "_participation_coefficient"), paste0("modularity_", suffix, "_within_module_degree"))
  mod_clu_prop$ID<-seq.int(nrow(mod_clu_prop))
  total_mod <- Reduce(function(...) merge(..., by='ID', all=T), list(mod_clu, mod_clu_role, mod_clu_prop))
  
  return(total_mod)
}

total_mod <- data.frame(node = '')
total_mod_phylo <- data.frame(node = '')
if (mod == T) {
  if (conventional == T) {
    total_mod <- Mod_analysis(paste0(taxon_dir, "output_mod/"), taxon_group, 'conv', "_tips_only_links_list_bipartite_network_")
  } 
  if (phylo == T && mod_phylo == T) {
    total_mod_phylo <- Mod_analysis(paste0(taxon_dir, "output_mod_phylo/"), taxon_group, 'phylo', "_all_branch_links_list_bipartite_network_")
  }
}

##########################################################
# load the geolocations
##########################################################
# make these myself...
# if (phylo == F) {
#   setwd(paste0(taxon_dir, "for_network_analyses_tips_only/"))
# } else {
#   setwd(paste0(taxon_dir, "for_network_analyses/"))
# }
# 
# locations <- paste0(taxon_group, "_network_geolayout.csv")
# locations <- read.csv(locations, header=F, skip=1)
# colnames(locations) <- c("ID.y", "x_coord",  "y_coord")
# # View(locations)


##########################################################
# load the data from biodiverse
##########################################################
setwd(paste0(taxon_dir))
biodiverse_export<- NA
# need to modify this bit.
# biodiverse_s2_cluster <- paste0(taxon_group, "_tablegroup_S2.csv")
if (do_biodiverse == T) {
  biodiverse_s2_cluster <- data.frame()
  biodiverse_Phylos2_cluster <- data.frame()
  if (conventional == T) {
    biodiverse_s2_cluster <- paste0(taxon_group, "_S2_tablegrouped.csv")
    # biodiverse_s2_cluster <- paste0(taxon_group, "_tablegroup_S2.csv")
    biodiverse_s2_cluster <- read.csv(biodiverse_s2_cluster, header=T)
  }
  if (phylo == T) {
    biodiverse_Phylos2_cluster <- paste0(taxon_group, "_PhyloS2_tablegrouped.csv")
    biodiverse_Phylos2_cluster <- read.csv(biodiverse_Phylos2_cluster, header=T)
  }
  
  # load the spatial analysis results
  biodiverse_spatial <- paste0(taxon_group, "_spatial.csv")
  # biodiverse_spatial <- paste0(taxon_group, "_spatial_biodiverse.csv")
  biodiverse_spatial <- read.csv(biodiverse_spatial, header=T)
  
  # biodiverse <- merge(biodiverse_s2_cluster, biodiverse_spatial)
  biodiverse_spatial$ELEMENT <- paste0("p_", biodiverse_spatial$ELEMENT)
  biodiverse_spatial$ELEMENT <- gsub(x = biodiverse_spatial$ELEMENT, pattern = "\\:", replacement = "_")
  
  biodiverse_export<-biodiverse_spatial$ELEMENT
  
  if (conventional == T) {
    biodiverse_export <- cbind(biodiverse_export, biodiverse_clusterS2 = as.numeric(biodiverse_s2_cluster$NAME))
  }
  if (phylo == T) {
    biodiverse_export <- cbind(biodiverse_export, biodiverse_phylo_clusterS2 = as.numeric(biodiverse_Phylos2_cluster$NAME))
  }
  
  biodiverse_export <- cbind(biodiverse_export, biodiverse_CWE = biodiverse_spatial$ENDW_CWE)
  biodiverse_export <- cbind(biodiverse_export, biodiverse_WE = biodiverse_spatial$ENDW_WE)
  biodiverse_export <- cbind(biodiverse_export, biodiverse_Richness = biodiverse_spatial$ENDW_RICHNESS)
  
  if (phylo == T){
    biodiverse_export <- cbind(biodiverse_export, biodiverse_PD = biodiverse_spatial$PD)
    biodiverse_export <- cbind(biodiverse_export, biodiverse_PE_CWE = biodiverse_spatial$PE_CWE)
    biodiverse_export <- cbind(biodiverse_export, biodiverse_PE_WE = biodiverse_spatial$PE_WE)
    biodiverse_export <- cbind(biodiverse_export, biodiverse_PHYLO_RARITY_CWR = biodiverse_spatial$PHYLO_RARITY_CWR)
    biodiverse_export <- cbind(biodiverse_export, biodiverse_RAREW_CWE = biodiverse_spatial$ENDW_RAREW_CWE)
    biodiverse_export <- cbind(biodiverse_export, biodiverse_RAREW_WE = biodiverse_spatial$ENDW_RAREW_WE)
  }
  
  biodiverse_export <- cbind(biodiverse_export, x_coord = biodiverse_spatial$Axis_0)
  biodiverse_export <- cbind(biodiverse_export, y_coord = biodiverse_spatial$Axis_1)
  
  colnames(biodiverse_export)[1] <- "node"
}
# View(locations)

##########################################################
# Calculating other network metrics
##########################################################
# need to create a function for these ones. Make it easier to change....


require(network)
library(igraph)
Network_analysis <- function(wd,path, suffix) {
  
  setwd(wd)
  myNetwork <- read.paj(path)
  
  
  vertexNames <- network.vertex.names(myNetwork)
  
  network_properties <- data.frame(node=vertexNames)
  
  myGraph_paj <- read_graph(path, format = "pajek")
  # plot(myGraph_paj)
  
  # network_properties$ID.y <- rownames(network_properties)
  prefix = paste0("network_", suffix, '_')
  network_properties[[paste0(prefix, "closeness")]] <- closeness(myGraph_paj)
  network_properties[[paste0(prefix, "betweenness")]] <- betweenness(myGraph_paj, directed = F)
  
  
  
  network_properties[[paste0(prefix, "degree")]] <- degree(myGraph_paj)
  network_properties[[paste0(prefix, "alpha_centrality")]] <- alpha_centrality(myGraph_paj)
  eigen_centrality <- eigen_centrality(myGraph_paj, directed = FALSE, scale = TRUE, weights = NULL, options = arpack_defaults)
  
  network_properties[[paste0(prefix, "eigen_centrality")]] <- eigen_centrality$vector
  
  return(network_properties)
}

network_properties <- data.frame(node = '')
network_phylo_properties <- data.frame(node = '')
network_rw_properties <- data.frame(node = '')
network_prw_properties <- data.frame(node = '')
if (do_networks == T) {
  if (conventional == T) {
    network_properties <- Network_analysis(paste0(taxon_dir, "for_network_analyses_tips_only/"), paste0(taxon_group, "_tips_only_pajek_no_weights.net"), 'conv')
  } 
  
  if (phylo == T) {
    network_phylo_properties <- Network_analysis(paste0(taxon_dir, "for_network_analyses/"), paste0(taxon_group, "_all_branch_pajek_single_branch_length.net"), 'phylo')
  }
  
  if (range_weighted == T) {
    network_rw_properties <- Network_analysis(paste0(taxon_dir, "for_network_analyses_tips_only/"), paste0(taxon_group, "_tips_only_pajek_range_weighted.net"), 'rw')
  }
  
  if (range_weighted_phylo == T) {
    network_prw_properties <- Network_analysis(paste0(taxon_dir, "for_network_analyses/"), paste0(taxon_group, "_all_branch_pajek_single_branch_length_range_weighted.net"), 'prw')
  }
}
network_all_properties <- Reduce(function(...) merge(..., all=T, by="node"), list(network_properties, network_phylo_properties, network_rw_properties, network_prw_properties))

# View(network_phylo_properties)
##########################################################
# putting it all together
##########################################################
empty <- data.frame(node = '')
for (i in 1:2) {
  if (i == 1) {
    locations_only = T
  } else {
    locations_only = F
  }
  total_map_all <- Reduce(function(...) merge(..., by='node', all=T), list(total_map, total_map_phylo, total_map_rw, total_map_prw))
  
  total_mod_all <- Reduce(function(...) merge(..., by='node', all=T), list( total_mod, total_mod_phylo))
  total_mod_all$ID.x <- NULL
  total_mod_all$ID.y <- NULL
  #View(total_map)
  
  total_map_all$trimnode <- strtrim(total_map_all$node, 23)
  total_mod_all$trimnode <- strtrim(total_mod_all$node, 23)
  
  total <- merge(total_map_all, total_mod_all, by="trimnode", all=T)
  
  total$node.y <- NULL
  names(total)[names(total) == 'node.x'] <- 'node'
  
  if (locations_only == T) {
    # total <- merge(total, locations, by="ID.y")
    if (do_networks == T) {
      total <- merge(total, network_all_properties, by="node")
    }
    total <- merge(total, biodiverse_export, by="node")
  } else {
    # total <- merge(total, locations, by="ID.y", all = TRUE)
    if (do_networks == T) {
      total <- merge(total, network_all_properties, by="node", all = TRUE)
    }
    total <- merge(total, biodiverse_export, by="node", all = TRUE)
  }
  
  
  # total$ID.x <- NULL
  # total$ID.y <- NULL
  # total$ID<-seq.int(nrow(total))
  # View(total)
  if (any(total$trimnode == '')) {
    total[total$trimnode == '',]$trimnode <- NA
  }
  total <- total[rowSums(is.na(total)) != ncol(total),]
  # total <- total[!duplicated(lapply(total, summary))]
  ##########################################
  # output 
  ##########################################
  setwd(taxon_dir)
  if (locations_only == T) {
  #   if (phylo == T) {
  #     write.csv(total, file=paste0(taxon_group,"_almalgamated_data_phylo_geo_only.csv"), row.names=FALSE)
  #     
  #   } else {
  #     write.csv(total, file=paste0(taxon_group,"_almalgamated_data_geo_only.csv"), row.names=FALSE)
  #     
  #   }
    write.csv(total, file=paste0(taxon_group,"_almalgamated_data_geo_only.csv"), row.names=FALSE)
  } else {
  #   if (phylo == T) {
  #     write.csv(total, file=paste0(taxon_group,"_almalgamated_data_phylo.csv"), row.names=FALSE)
  #     
  #   } else {
  #     write.csv(total, file=paste0(taxon_group,"_almalgamated_data.csv"), row.names=FALSE)
  #     
  #   }
    write.csv(total, file=paste0(taxon_group,"_almalgamated_data.csv"), row.names=FALSE)
  }
}
