#################################################
# This file loads a tree file and biodiverse label export then trims them 
# and generates the files needed to create a network analysis

# to use this, you need to load the tree into figtree, and export it. Otherwise they don't like the formatting.
#################################################
library(ape)
library(phytools)
library(tools)

taxon_group <- "fernsPhylo" 
output_folder <- paste0("Z:/!my_working/data/", taxon_group,  "/for_network_analyses/")
dir.create(output_folder, showWarnings = FALSE)
setwd(paste0("Z:/!my_working/data/", taxon_group, "/"))

# export tree from biodiverse then save from figtree. Make sure is newick export.
tree_file  <- paste0(taxon_group, "_clean_tree_figtree.nwk")
tree  <- read.newick(file=tree_file)
tree  <- collapse.singles(tree, root.edge = FALSE) # need to do this because of Biodiverse pruning I think

biodiverse_label_export_file <- paste0(taxon_group, "_label_export.csv")

##########################################
#visualise loaded tree
# plot(tree, root.edge=TRUE)
# this lines up tip numbers & labels for simplicity
# par(mar=c(1,1,1,1)) # set margins
# nodelabels(adj = c(0.2, 0.2), font = 1, cex=.6,  bg = "red")
# edgelabels(font = 1, cex=.6, bg = "white")#adj = c(-1, -1),  
#is.rooted(tree)

#####################################
#load tree and prune to orrurance matrix data loaded from biodiverse label export
#####################################
#basename(tree_file)#filename only
#file_path_sans_ext(tree_file) #file name only with path
#basename(file_path_sans_ext(tree_file))#file name only with no extension
print("loading tree and occurrence matix and pruning")

output_occurrence_matrix_file <-   paste(output_folder, taxon_group, "_", basename(file_path_sans_ext(biodiverse_label_export_file)),"_trimmed_to_tree", ".csv",sep="")
output_tree_file <-  paste(output_folder, taxon_group, "_", basename(file_path_sans_ext(tree_file)),"_pruned_to_spatial_data", ".new",sep="")# auto replace last 4 characters with new extension and add projection to and of filename

pt <- read.csv(file=biodiverse_label_export_file, header=T,sep=",",check.names=FALSE)# load file
# View(pt)
#ncol(pt)
occurrence_matrix_df <- as.data.frame(pt) #all spatial data
trimmed_occurrence_matrix <- occurrence_matrix_df[occurrence_matrix_df$Axis_0 %in% tree$tip.label, ]
# View(trimmed_occurrence_matrix)
#data not in matrix, but in tree
deleted_taxa_occurrence_matrix <- occurrence_matrix_df[ ! occurrence_matrix_df$Axis_0 %in% tree$tip.label, ]
# View(deleted_taxa_occurrence_matrix)

all_tree_labels <- as.data.frame(tree$tip.label)
colnames(all_tree_labels) <- c("taxa")
#View(all_tree_labels)
in_tree_but_not_in_matrix <- subset(all_tree_labels, !(all_tree_labels$taxa  %in% occurrence_matrix_df$Axis_0))

# View(in_tree_but_not_in_matrix)
tips_to_drop <- as.vector(in_tree_but_not_in_matrix$taxa)
pruned_tree <- drop.tip(tree, tips_to_drop)
#pruned_tree$tip.label
plot(pruned_tree)

write.tree(pruned_tree, file = output_tree_file, append = FALSE, digits = 10, tree.names = FALSE)

################################################################################
#View(pt)
#ncol(pt)

colnames(trimmed_occurrence_matrix) <- paste("p", colnames(trimmed_occurrence_matrix), sep="_")
#names(trimmed_occurance_matrix) <- gsub(x = names(trimmed_occurance_matrix), pattern = "\\:", replacement = "|")
names(trimmed_occurrence_matrix) <- gsub(x = names(trimmed_occurrence_matrix), pattern = "\\:", replacement = "_")
names(trimmed_occurrence_matrix)[2]<-"taxa"

locs_only <- trimmed_occurrence_matrix[3:ncol(trimmed_occurrence_matrix)] 
#View(locs_only)
locs_only[is.na(locs_only)] <- 0 # replace the NA's with 0's
#View(taxa_only)
locs_only[locs_only > 0.1] <- 1 # replace all value greater than 0.1 with 1
#View(locs_only)

occurrence_matrix_taxa <- trimmed_occurrence_matrix[2] # get the axis columns to reproject back to lat long

###################################################
#reproject file if needed (needs to be group export from biodiverse not label export)
#coordinates(pt) <- c("Axis_0",  "Axis_1")    # set spatial coordinates
#wgs84proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #WGS84 proj4 projection
#albersGDA94proj <- "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" #GDA94 / Australian Albers proj4 projection
#projection(pt) <- CRS(albersGDA94proj)
#pt <- spTransform(pt, CRS(wgs84proj))
#projection(pt)
#View(pt)
################################################
final_occurrence_matrix <- cbind(occurrence_matrix_taxa, locs_only)
write.csv(final_occurrence_matrix, file=output_occurrence_matrix_file, row.names = FALSE)

# View(final_occurrence_matrix)
#################################################

tree <- pruned_tree
root_number <- length(tree$tip.label)+1 #root is N taxa + 1
print(paste0("pruned tree contains: ", length(tree$tip.label), " tips"))
print(paste0("pruned tree written: ", output_tree_file))
print(paste0("trimmed occurrence matrix: ", output_occurrence_matrix_file))

###################################################

master_branch_list <- data.frame(tree$edge)
master_branch_list$branch_length <- tree$edge.length
master_branch_list$branch_number <- rownames(master_branch_list)
master_branch_list$branch_name <- paste("b_00", rownames(master_branch_list), sep="")

names(master_branch_list)[1]<-"left_node"
names(master_branch_list)[2]<-"right_node"

tips <- data.frame(branch_name = tree$tip.label)
tips$right_node <- rownames(tips)

master_branch_list <- merge(x = master_branch_list, y = tips, by.x = "right_node", by.y="right_node", all.x=TRUE)

ColsToMerge <- c("branch_name.y", "branch_name.x")
master_branch_list[["merged_branch_names"]] <-  apply(master_branch_list[, ColsToMerge], 1, function(rr) ifelse(is.na(rr[[1]]), rr[[2]], rr[[1]]))# merge the 2 coulmns with branch names so tips take precedence then branch number


###################
# Loop through the taxa in the tree and generate a branch_by_taxa file tracing each taxa back to the root indicating the branch numbers
####################
taxa_by_branches <- data.frame(taxa=tree$tip.label) 

for(i in 1:length(tree$tip.label)){
  tip_to_trace <- i
  nodes_traced <- nodepath(tree, from = root_number, to = tip_to_trace)
  branches <- subset(master_branch_list, master_branch_list[,"right_node"] %in% nodes_traced & master_branch_list[,"left_node"] %in% nodes_traced)
  branches_numbers_as_list <- list(branches$branch_number)
  # print(branches_numbers_as_list)
  taxa_by_branches$branches[i] <- list(unlist(branches_numbers_as_list))
}
taxa_by_branches$branches <- vapply(taxa_by_branches$branches, paste, collapse = ", ", character(1L))

#View(taxa_by_branches)
#setwd("C:/GIS-Datasets/Phylogenetic_network")
#write.csv(taxa_by_branches, file=paste(output_folder, taxon_group, "_", "by_branches.csv", sep=""), row.names = FALSE)#write taxa list with each branch tracing to root

########################################
#loop through the master branch list and calculate a table with the summed branch lengths from the root to the branch for all branches including the tips
#######################################
All_branches_with_summed_lengths <- data.frame(taxa=master_branch_list$merged_branch_names) 
for(i in 1:nrow(master_branch_list)){
     #print(master_branch_list[i,]$right_node)
     node_to_trace <- master_branch_list[i,]$right_node
     nodes_traced <- nodepath(tree, from = root_number, to = node_to_trace)
     branches <- subset(master_branch_list, master_branch_list[,"right_node"] %in% nodes_traced & master_branch_list[,"left_node"] %in% nodes_traced)
     branches_numbers_as_list <- list(branches$branch_number)
     #print(branches_numbers_as_list)
    All_branches_with_summed_lengths$branches[i] <- list(unlist(branches_numbers_as_list))
  for(p in All_branches_with_summed_lengths$branches[i]){
      summed_weight <- 0
       for(k in p){
         # print(k)
         branch_row <- master_branch_list[ which(master_branch_list$branch_number==k),] 
         summed_weight <- summed_weight + branch_row$branch_length
#          print(branch_row$branch_length)
#          print(summed_weight)
       }
       All_branches_with_summed_lengths$summed_branch_lengths[i] <- summed_weight
       }
}

All_branches_with_summed_lengths$branches <- vapply(All_branches_with_summed_lengths$branches, paste, collapse = ", ", character(1L))

colnames(All_branches_with_summed_lengths) <- c("branch_name_for_checking", "branch_numbers_to_root", "summed_branch_length") 
#View(All_branches_with_summed_lengths)
#write.csv(All_branches_with_summed_lengths, file=paste(output_folder, taxon_group, "_", "all_branches_tracing_branch_lengths.csv", sep=""), row.names = FALSE)#write taxa list with each branch tracing to root

master_branch_list <- cbind(master_branch_list, All_branches_with_summed_lengths)

tab_tree_filename <-  paste(output_folder,basename(file_path_sans_ext(output_tree_file)), "_tabular_tree.csv", sep="")

write.csv(master_branch_list, file=tab_tree_filename, row.names = FALSE)#write taxa list with each branch tracing to root
#View(master_branch_list)
print(paste0("Tabular tree written: ", tab_tree_filename))


########################

############
#calculate the subtending tips to any given branch in the tree and add it to the master tabular tree  
#

master_branch_list$subtending_tips <- 0

for(i in 1:nrow(master_branch_list)) {
  node_to_trace <- master_branch_list[i,"right_node"] #loop through the right node column
  descendant_nodes <- getDescendants(tree, node_to_trace, curr=NULL) # get the decents
  # print(descendant_nodes)
  master_branch_list[i,"subtending_tips"] <- length(which(descendant_nodes <= length(tree$tip.label))) # count all the ones that are tips and add it to a new column
}

# View(master_branch_list)

write.csv(master_branch_list, file=tab_tree_filename, row.names = FALSE)#write taxa list with each branch tracing to root
#View(master_branch_list)
#############################################################









####################################
#calc evolutionary distinctiveness NOT WORKING WITH OTHER PACKAGES ETC. 
# #
# library(picante)
# #evol.distinct(tree, type = c("equal.splits", "fair.proportion"), scale = FALSE, use.branch.lengths = TRUE)
# evol_dist <- evol.distinct(tree, type = "equal.splits", scale = FALSE, use.branch.lengths = TRUE)
# ?evol.distinct
# 
# 
# biodiverse_label_export_file <- "biodiverse_label_export.csv"
# pt <- read.csv(file=biodiverse_label_export_file, header=T,sep=",",check.names=FALSE)# load file
# row.names(pt) <- pt[,2] 
# pt <- pt[3:ncol(pt)]
# #View(locs_only)
# pt[is.na(pt)] <- 0 # replace the NA's with 0's
# #View(taxa_only)
# pt[pt > 0.1] <- 1 # replace all value greater than 0.1 with 1
# #match.phylo.data(phy, data)
# test_picante_prune <- match.phylo.data(tree, pt) 
# # View(pt)
# str(test_picante_prune)
# 
# 
# #
# #test<- pd("cyrtoloba", tree, include.root=TRUE)
# 
# #plot(tree)
# #tree$tip.label
# #str(tree)
# 
# # library(devtools)
# # install_git("https://github.com/eliotmiller/ecoPDcorr.git")
# # library(ecoPDcorr)
# # test_mat <- "C:/GIS-Datasets/Phylogenetic_network/test_ed_weeds.csv"
# # pt <- read.csv(file=test_mat, header=T,sep=",",check.names=FALSE)# load file
# # 
# # View(pt)
# # 
# # rownames(pt) <- pt[,c(1)]
# # pt <- pt[, -c(1)] #drop unwanted columns
# # #pt <- as.matrix(pt)
# # 
# # plot(weeds.tree)
# # str(weeds.tree)
# # 
# # ?read.csv
# # 
# # test_phylo <- phylo4com(weeds.tree, pt)
# # 
# # 
# # aed(test_phylo)
# # ed(test_phylo)
# # haed(test_phylo)
# # eaed(test_phylo)
# # 
# 
# 
# ####################
# #calculate range of each taxa in final matrix
# # #####
# # library(dplyr)
# # final_occurrence_dt <- tbl_df(final_occurrence_matrix) #turn the matrix into a dplyr data table
# # taxon_ranges <- final_occurrence_dt
# # 
# # taxon_ranges <- final_occurrence_dt %>%
#   mutate(range = rowSums(.[2:ncol(taxon_ranges)])) %>% #add a range column that sums the row
#   select(taxa, range) # just select the columns we want
# 
# ranges_occurrence_matrix_file <-   paste(output_folder, taxon_group, "_", basename(file_path_sans_ext(biodiverse_label_export_file)),"_trimmed_to_tree_ranges", ".csv",sep="")
# #write.csv(taxon_ranges, file=ranges_occurrence_matrix_file, row.names = FALSE)
# View(taxon_ranges)
# master_branch_list <- merge(x = master_branch_list, y = taxon_ranges, by.x = "branch_name.y", by.y="taxa", all.x=TRUE)#add it to the tabular tree
# View(master_branch_list)

