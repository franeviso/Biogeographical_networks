#################################################################################
# This file takes the occurance_matrix and tabular tree and generates all the 
# network file for running network anaylses.
# pajek files:
#             links_list_bipartite_network.txt            (modularity analysis)
#             pajek_single_branch_length.net              (map equation analysis)
#             pajek_summed_branch_length.net              (map equation analysis)
#             tips_only_links_list_bipartite_network.txt  (map equation analysis)
#             conventional_pajek_no_weights.net           (map equation analysis)
#################################################################################
library(bipartite)

# rscript "Y:\nathaniel_load_files_and_create_network_matrix.R" AustralianGeneraPhylo
 
#input files
# range_weighted <- T
args <- commandArgs(trailingOnly = TRUE) 
taxon_group <- args[1] 
# taxon_group <- "fernsPhylo" 
output_folder <- paste0("Z:/!my_working/data/", taxon_group,  "/for_network_analyses/")
dir.create(output_folder, showWarnings = FALSE)

setwd(output_folder)
taxa_tabular_tree_export_file  <- paste0(taxon_group,"_",taxon_group,"_clean_tree_figtree_pruned_to_spatial_data_tabular_tree.csv")
co_occurance_file <- paste0(taxon_group,"_",taxon_group,"_label_export_trimmed_to_tree.csv")
  # "eucalypts_2015Phylo_eucalypts_2015Phylo_label_export_trimmed_to_tree.csv"
# taxa_tabular_tree_export_file  <- "glycine_glycine_new_exported_from_biodiverse_pruned_to_spatial_data_tabular_tree.csv"
# co_occurance_file <- "glycine_glycine_trimmed_to_tree.csv"

######################################

taxa_tabular_tree_export  <- read.csv(file= taxa_tabular_tree_export_file, header=T,sep=",", check.names=FALSE)
#View(taxa_tabular_tree_export)
#nrow(taxa_tabular_tree_export)

branch_name <- data.frame(branch_number=taxa_tabular_tree_export$branch_number,  branch_name=taxa_tabular_tree_export$merged_branch_names) # create branch_name list to use down  below 
branch_name <-  branch_name[with(branch_name, order(branch_number)), ]
#View(branch_name)

taxa_list_matrix  <- taxa_tabular_tree_export[complete.cases(taxa_tabular_tree_export),]# trim off internal branches
#View(taxa_list_matrix)
taxa_list_matrix_asc <- taxa_list_matrix[order(taxa_list_matrix$merged_branch_names),]
taxa_list_matrix_asc <- as.matrix(taxa_list_matrix_asc)
#View(taxa_list_matrix_asc)
#str(taxa_list_matrix_asc)

#load co-occurance matrix #group export from biodiverse
co_occurance_matrix  <- read.csv(file= co_occurance_file, header=T,sep=",", check.names=FALSE)
#View(co_occurance_matrix)
#ncol(co_occurance_matrix)
#colnames(co_occurance_matrix)

locations_only  <- colnames(co_occurance_matrix[,-1])
#locations_only

#create blank matrix
network_matrix <- matrix (0,nrow=nrow(taxa_tabular_tree_export), ncol=ncol(co_occurance_matrix))
colnames(network_matrix) <- colnames(co_occurance_matrix)
colnames(network_matrix)[1] <- "branch"
network_matrix[,1] <- seq(1:nrow(taxa_tabular_tree_export))#create a unique number for each branch and add to first column
#View(network_matrix)
#dim(network_matrix)

print("Making network matrix")
d <- 1 # counter for locations in the columns in matrix
for(i in co_occurance_matrix){
 #print(i)
   c <- 1 # counter for taxa rows in taxa list
  for(j in i){
      #print(j)
    if(j == 1){
      #print(taxa_list_matrix_asc[c,2])
      #branch_list <- as.list(taxa_list_matrix_asc[c,2]) #get the list of branches a taxa contains
      branch_list <- as.list(taxa_list_matrix_asc[c,"branch_numbers_to_root"]) #get the list of branches a taxa contains
      for(k in branch_list){
        branch <- as.vector(unlist(strsplit(k,",")),mode="list") # split it into single branches
        for(p in branch){
         # print(p)
         # print(i)
          network_matrix[as.numeric(p),d] <-1 #populate the newtork matrix cells with 1 for each branch location 
        }
      }
     }   
    c <- c+1
  }
  d <- d+1
}
#View(network_matrix)

write.csv(network_matrix, file=paste(output_folder, taxon_group, "_all_branch_occurance_matrix.csv", sep=""), row.names = FALSE)#write newtok matrix

network_matrix_only<- network_matrix[,-1] 
#View(network_matrix_only)

#locations_only_file  <- "C:/GIS-Datasets/Glycine/locations_only.csv"
#locations_only  <- read.csv(file=locations_only_file, header=T,sep=","check.names=FALSE)
#View(locations_only)
locs <- unlist(locations_only)
#locs
#locations_only$locations[4]

print("making links list network file")

dim(network_matrix_only)[1]
con_list <-c(0,0)
for (i in 1:dim(network_matrix_only)[1]){# i is branches
  #print(i)
  for (j in 1:dim(network_matrix_only)[2]){
    if(network_matrix_only[i,j] == 1){ # j is locations
      #con_list <- rbind(con_list, c((dim(network_matrix_only)[2]+i),j)) # original number only version
      con_list <- rbind(con_list, c(paste0(branch_name$branch_name[i]),paste0(locs[j])))#substitute strings with locations and sp names
    }
  }
}
#branch_name[2]
#View(con_list)
con_list <- con_list[-1,]


output_network <- con_list
output_network[,1] <- strtrim(con_list[,1], 23)
write.table(output_network, file=paste(output_folder, taxon_group, "_all_branch_links_list_bipartite_network.txt", sep=""), row.names = FALSE, col.names=FALSE,  quote = FALSE, sep=" ")#write newtork matrix, space delimited text file used for modularity analyses


#View(con_list)
#write.csv(con_list, file=paste(output_folder, taxon_group, "_links_list_bipartite_network.csv", sep=""), row.names = FALSE)#write newtork matrix
# write.table(con_list, file=paste(output_folder, taxon_group, "_all_branch_links_list_bipartite_network.txt", sep=""), row.names = FALSE, col.names=FALSE,  quote = FALSE, sep=" ")#write newtork matrix, space delimited text file used for modularity analyses

###############################################################
# write pajek format
###############################################################
total_nodes <- nrow(branch_name) + length(locations_only)
verticies_line <- paste0("*Vertices ", total_nodes)
########
verticies_geo <- data.frame(vert_name=locations_only)
verticies_phylo <- data.frame(vert_name=branch_name$branch_name)
verticies_all <- rbind(verticies_geo, verticies_phylo)
verticies_all <- cbind(vert_number = rownames(verticies_all), verticies_all) # add the rownames as a proper column
#View(verticies_all)

verticies_geo_with_numbers <- verticies_all[c(1:nrow(verticies_geo)),] 
#View(verticies_geo_with_numbers)

verticies_phylo_with_numbers <- verticies_all[-(1:nrow(verticies_geo)),] 
#View(verticies_phylo_with_numbers)

#order the tree table by branch number so that the correct summed branch lengths can be added
taxa_tabular_tree_export_odered_by_branch_number <- taxa_tabular_tree_export[with(taxa_tabular_tree_export, order(branch_number)), ]
#View(taxa_tabular_tree_export_odered_by_branch_number)

#############################################################
# make pajek edge list with single branch lengths as weight
#############################################################
print("making pajek edge list with single branch lengths as weight")
dim(network_matrix_only)[1]
edge_list <- data.frame(edge1= integer(0), edge2= integer(0), weight = character(0))
#edge_list <- data.frame()
for (i in 1:dim(network_matrix_only)[1]){# i is branches
  #print(i)
  for (j in 1:dim(network_matrix_only)[2]){
    if(network_matrix_only[i,j] == 1){ # j is locations
      #con_list <- rbind(con_list, c(i,j)) # original number only version
      #print(as.character(verticies_phylo_with_numbers$vert_number[i]))
      #print(j)
      new_data <- data.frame(edge1= verticies_phylo_with_numbers$vert_number[i], edge2= verticies_geo_with_numbers$vert_number[j], weight =taxa_tabular_tree_export_odered_by_branch_number$branch_length[i])# generates individual branch length as weight
      #new_data <- data.frame(edge1= verticies_phylo_with_numbers$vert_number[i], edge2= verticies_geo_with_numbers$vert_number[j], weight =taxa_tabular_tree_export$summed_branch_length[i])# generates summed branch length as weight
      edge_list <- rbind(edge_list, new_data)#substitute strings with locations and sp names
    }
  }
}
#View(edge_list)
#write.csv(edge_list, file="pajek_links_list.csv", row.names = FALSE)#write newtok matrix
##########################################

total_links <- nrow(con_list)
edges <- paste0("*Edges ", total_links) 

verticies_phylo_range <- data.frame(rowSums(network_matrix_only))
verticies_phylo_range$ID <- verticies_phylo_with_numbers$vert_number

for(i in 1:2) {
  if (i == 2) {
    filename <- paste(output_folder, taxon_group, "_all_branch_pajek_single_branch_length_range_weighted.net", sep="")
    old_edge_list <- edge_list
#     for(j in 1:nrow(edge_list)){
#       edge_list$weight[j] <- edge_list$weight[j]/(nrow(edge_list[edge_list[, "edge1"] == edge_list$edge1[j],]))
#     }
    edge_list$weight <- edge_list$weight/verticies_phylo_range[as.numeric(as.character(edge_list$edge1))- min(as.numeric(as.character(edge_list$edge1))) + 1, 1]
  } else {
    filename <- paste(output_folder, taxon_group, "_all_branch_pajek_single_branch_length.net", sep="")
    
  }
  #cat("# Pajek format")
  # cat("\n")
  sink(filename)
  cat(verticies_line)
  cat("\n")
  for(i in 1:nrow(verticies_all)){
    cat(paste(i, " ", "\"", verticies_all$vert_name[i], "\" ", sep=""))
    cat("\n")
  }
  cat(edges)
  cat("\n")
  sink()
  # cat("\n")
  # for(i in 1:nrow(edge_list)){
  #   if (range_weighted == T) {
  #     cat(paste(edge_list$edge1[i], edge_list$edge2[i], edge_list$weight[i]/(nrow(edge_list[edge_list[, "edge1"] == edge_list$edge1[i],])), sep=" "))
  #   } else {
  #     cat(paste(edge_list$edge1[i], edge_list$edge2[i], edge_list$weight[i], sep=" "))
  #   }
  #   cat("\n")
  # }
  # sink()
  write.table(edge_list, file=filename, col.names = F, row.names = F, quote=F, append=T)
}
##########################################