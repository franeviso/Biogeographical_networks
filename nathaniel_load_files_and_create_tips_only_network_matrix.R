#################################################################################
# This file takes the occurance_matrix only and generates all the 
# network file for running network anaylses.
# pajek files:
#             tips_only_links_list_bipartite_network.txt  (modularity)
#             conventional_pajek_no_weights.net   (map equation analysis)
#             _network_geolayout.csv              (for gephi geolayout)
#             _network_branch_label_only.csv      (for labeling only tips/branches)
#             _network_properties_tips_only.csv   (calc network properties to use for correlations)
#
#################################################################################
library(bipartite)
library(tools)

#input files
taxon_group <- "fernsPhylo" 
range_weight <- F
output_folder <- paste0("Z:/!my_working/data/", taxon_group,  "/for_network_analyses_tips_only/")
dir.create(output_folder, showWarnings = FALSE)
#setwd(output_folder)
setwd(paste0("Z:/!my_working/data/", taxon_group, "/"))
biodiverse_label_export_file <- paste0(taxon_group, "_label_export.csv")

output_occurrence_matrix_file <-   paste(output_folder, taxon_group, "_", basename(file_path_sans_ext(biodiverse_label_export_file)),"_occurrence_matrix_tips_only_names", ".csv",sep="")

pt <- read.csv(file=biodiverse_label_export_file, header=T,sep=",",check.names=FALSE)# load file
# View(pt)
#ncol(pt)
occurrence_matrix_df <- as.data.frame(pt) #all spatial data

################################################################################
#View(pt)
#ncol(pt)

colnames(occurrence_matrix_df) <- paste("p", colnames(occurrence_matrix_df), sep="_") #add a p_ to the front of colnames
names(occurrence_matrix_df) <- gsub(x = names(occurrence_matrix_df), pattern = "\\:", replacement = "_") 
names(occurrence_matrix_df)[2]<-"taxa"
#View(occurrence_matrix_df)
locs_only <- occurrence_matrix_df[3:ncol(occurrence_matrix_df)] 
#View(locs_only)
locs_only[is.na(locs_only)] <- 0 # replace the NA's with 0's
#View(taxa_only)
locs_only[locs_only > 0.1] <- 1 # replace all value greater than 0.1 with 1
#View(locs_only)

occurrence_matrix_taxa <- occurrence_matrix_df[2] # get the axis columns to reproject back to lat long

final_occurrence_matrix <- cbind(occurrence_matrix_taxa, locs_only)

# View(final_occurrence_matrix)
write.csv(final_occurrence_matrix, file=output_occurrence_matrix_file, row.names = FALSE)

##########################################################################################
# Generate conventional tip only pajek files for comparative analysis
##########################################################################################

conventional_branch_name <- data.frame(branch_name=final_occurrence_matrix$taxa) # create branch_name l
#conventional_branch_name <- as.data.frame(conventional_branch_name[order(conventional_branch_name$branch_name, decreasing=FALSE), ])
conventional_branch_name <- cbind(branch_number = rownames(conventional_branch_name), conventional_branch_name) # add the rownames as a proper column

colnames(conventional_branch_name) <- c("branch_number", "branch_name")
# View(conventional_branch_name)

#View(co_occurance_matrix)

#create blank matrix
conventional_network_matrix <- final_occurrence_matrix
colnames(conventional_network_matrix)[1] <- "branch"
conventional_network_matrix[,1] <- seq(1:nrow(conventional_branch_name))#create a unique number for each branch and add to first column

#View(conventional_network_matrix)

write.csv(conventional_network_matrix, file=paste(output_folder, taxon_group, "_tips_only_occurance_matrix_numbers.csv", sep=""), row.names = FALSE)#write newtok matrix

conventional_network_matrix_only <- conventional_network_matrix[,-1] 
#View(conventional_network_matrix_only)

locs <- unlist(colnames(conventional_network_matrix_only))
locs

print("making links list conventional_network_matrix file")

dim(conventional_network_matrix_only)[1]
con_list <-c(0,0)
for (i in 1:dim(conventional_network_matrix_only)[1]){# i is branches
  #print(i)
  for (j in 1:dim(conventional_network_matrix_only)[2]){
    if(conventional_network_matrix_only[i,j] == 1){ # j is locations
      #con_list <- rbind(con_list, c((dim(conventional_network_matrix_only)[2]+i),j)) # original number only version
      con_list <- rbind(con_list, c(paste0(conventional_branch_name$branch_name[i]),paste0(locs[j])))#substitute strings with locations and sp names
    }
  }
}
# View(conventional_network_matrix_only)
# View(con_list)
con_list <- con_list[-1,]
#View(con_list)
#write.csv(con_list, file=paste(output_folder, taxon_group, "_links_list_bipartite_network.csv", sep=""), row.names = FALSE)#write newtork matrix

output_network <- con_list
output_network[,1] <- strtrim(con_list[,1], 23)
write.table(output_network, file=paste(output_folder, taxon_group, "_tips_only_links_list_bipartite_network.txt", sep=""), row.names = FALSE, col.names=FALSE,  quote = FALSE, sep=" ")#write newtork matrix, space delimited text file used for modularity analyses


###############################################################
# write pajek format conentional tips only
###############################################################
total_nodes <- nrow(conventional_branch_name) + length(locs)
verticies_line <- paste0("*Vertices ", total_nodes)
########
verticies_geo <- data.frame(vert_name=locs)
verticies_phylo <- data.frame(vert_name=conventional_branch_name$branch_name)
verticies_all <- rbind(verticies_geo, verticies_phylo)
verticies_all <- cbind(vert_number = rownames(verticies_all), verticies_all) # add the rownames as a proper column
#View(verticies_all)
#View(verticies_phylo_range)
verticies_geo_with_numbers <- verticies_all[c(1:nrow(verticies_geo)),] 
#View(verticies_geo_with_numbers)

verticies_phylo_with_numbers <- verticies_all[-(1:nrow(verticies_geo)),] 
#View(verticies_phylo_with_numbers)



verticies_phylo_range <- data.frame(rowSums(conventional_network_matrix_only))
verticies_phylo_range$ID <- verticies_phylo_with_numbers$vert_number
# View(verticies_phylo_range)
#order the tree table by branch number so that the correct summed branch lengths can be added
#taxa_tabular_tree_export_odered_by_branch_number <- taxa_tabular_tree_export[with(taxa_tabular_tree_export, order(branch_number)), ]
#View(taxa_tabular_tree_export_odered_by_branch_number)

#################################
#################################
# make pajek edge list with no weight
#################################
print("making conventional pajek edge list with no weight")
dim(conventional_network_matrix_only)[1]
edge_list <- data.frame(edge1= integer(0), edge2= integer(0))
#edge_list <- data.frame()
for (i in 1:dim(conventional_network_matrix_only)[1]){# i is branches
  #print(i)
  for (j in 1:dim(conventional_network_matrix_only)[2]){
    if(conventional_network_matrix_only[i,j] == 1){ # j is locations
      #con_list <- rbind(con_list, c(i,j)) # original number only version
      #print(as.character(verticies_phylo_with_numbers$vert_number[i]))
      #print(j)
      new_data <- data.frame(edge1= verticies_phylo_with_numbers$vert_number[i], edge2= verticies_geo_with_numbers$vert_number[j])# generates individual branch length as weight
      #new_data <- data.frame(edge1= verticies_phylo_with_numbers$vert_number[i], edge2= verticies_geo_with_numbers$vert_number[j], weight =taxa_tabular_tree_export$summed_branch_length[i])# generates summed branch length as weight
      edge_list <- rbind(edge_list, new_data)#substitute strings with locations and sp names
    }
  }
}
# edge_list$weight <- NULL
#write file
total_links <- nrow(con_list)
edges <- paste0("*Edges ", total_links) 
for (i in 1:2) {
  if (i == 2) {
    filename <- paste(output_folder, taxon_group, "_tips_only_pajek_range_weighted.net", sep="")
    old_edge_list <- edge_list
    edge_list$weight <- 1/verticies_phylo_range[as.numeric(as.character(edge_list$edge1))- min(as.numeric(as.character(edge_list$edge1))) + 1, 1]

  } else {
    filename <- paste(output_folder, taxon_group, "_tips_only_pajek_no_weights.net", sep="")
  }
  sink(filename)
  #cat("# Pajek format")
  #cat("\n")
  cat(verticies_line)
  cat("\n")
  for(i in 1:nrow(verticies_all)){
    cat(paste(i, " ", "\"", verticies_all$vert_name[i], "\" ", sep=""))
    cat("\n")
  }
  cat(edges)
  cat("\n")
#   for(i in 1:nrow(edge_list)){
#     if (range_weight == T) {
#       cat(paste(edge_list$edge1[i], edge_list$edge2[i], 1/verticies_phylo_range[as.numeric(as.character(edge_list$edge1[i]))- min(as.numeric(as.character(edge_list$edge1))) + 1, 1], sep=" "))
#     } else {
#       cat(paste(edge_list$edge1[i], edge_list$edge2[i],  sep=" "))
#     }
#     cat("\n")
#   }
  sink()
  
  write.table(edge_list, file=filename, col.names = F, row.names = F, quote=F, append=T)
  
}
# View(verticies_phylo_range)
##########################################


