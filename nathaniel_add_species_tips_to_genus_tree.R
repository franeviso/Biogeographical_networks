library(ape)
library(phytools)
library(tools)       


bind.tip<-function(tree,tip,where=NULL){
  obj<-bind.tree(tree,tip,where=where)
  return(obj)
}

taxon_group <- "fernsPhylo" 
setwd(paste0("Y:/data/", taxon_group, "/"))
tree_file  <- paste0("RAxML_bestTree.rerooted.tre")
species_file <- "updated_AVH_ferns_09_14_SPECIESLIST.csv"
tree  <- read.tree(file=tree_file)
par(mar=c(1,1,1,1)) # set margins
plot(tree)
nodelabels(adj = c(0.2, 0.2), font = 1, cex=.6,  bg = "red")
edgelabels(font = 1, cex=.6, bg = "white")#adj = c(-1, -1),

# str(tree)
# View(tree$tip)
# View(as.data.frame(tree$tip.label))
# View(tree$edge)
# create tip (modify as desired)

species_list <- read.table(file=species_file)
genus_species = as.data.frame(within(species_list, V1<-data.frame(do.call('rbind', strsplit(as.character(V1), '_', fixed=TRUE)))))
# View(genus_species)
species_list$genus <- as.data.frame(genus_species[, 1][1])
species_list$species <- as.data.frame(genus_species[, 1][2])

# colnames(species_list) <- c("genus_species", "genus", "species")
# View(species_list)


master_branch_list <- data.frame(tree$edge)
master_branch_list$branch_length <- tree$edge.length

length <- mean(master_branch_list$branch_length)

i <- 1
oldtree <- tree
for (genus in tree$tip.label) {
  print(i)
  print(genus)
  
  species_subset <- species_list[species_list[, "genus"] == genus, "species"]
  
  tree_string <- '('
  # m = 1
  for (species in species_subset[, 1]) {
    # print(m)
    # print(species)
    # print(tree_string)
    tree_string <- paste0(tree_string, genus, '_', species, ":", length, ",")
    # m = m + 1
  }
  
  tree_string <- paste0(substr(tree_string, 1, nchar(tree_string)-1), ');')
  tip_tree  <- read.tree(text=tree_string)
  # plot(tip_tree)
  tree <- bind.tree(tree,tip_tree,where=i)
  
  i = i + 1
}

plot(tree)
write.tree(tree, file = paste0(tree_file, "_species_tree.nwk"), append = FALSE, digits = 10, tree.names = FALSE)
##################



