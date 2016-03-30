
library(ape)
library(phytools)
library(tools)

taxon_group <- "acaciaPhylo"
taxon_dir <- paste0("Z:/!my_working/data/", taxon_group, "/")
dendro_type <- "S2"

tree_file <- paste0(taxon_dir, taxon_group, '_', dendro_type, '_tablegrouped_dendro.nwk')
output_folder <- paste0(taxon_dir, 'networks/')
dir.create(output_folder, showWarnings = FALSE)
# print(tree_file)

String <- paste(readLines(tree_file))
String <- gsub("'([^:_]+)\\:([^']+)'","\\1\\_\\2",String)
String <- substr(String, 1, nchar(String)-2)
String <- paste0(String, ';') # need to remove the :0 at the end of hte biodiverse newick export and add a ; instead.
# print(String)
new_file <- paste0(tree_file, '.changed')
# print(new_file)
write(String, file = new_file)


tree  <- read.tree(file=new_file)
# str(tree)
# summary(tree)
# plot(tree, root.edge=TRUE)


# color groupings...
group_file <- paste0(taxon_dir, taxon_group, '_almalgamated_data_geo_only.csv')

to_group <-read.csv(group_file, header=T,sep=",") # read data to generate raster files
if (dendro_type == "S2") {
  all_cells <- subset(to_group,  select=c(x_coord, y_coord, biodiverse_clusterS2)) # r generate a rater based on the COLUORS column
  all_cells$NAME <- as.numeric(all_cells$biodiverse_clusterS2)
} else {
  all_cells <- subset(to_group,  select=c(x_coord, y_coord, biodiverse_phylo_clusterS2)) # r generate a rater based on the COLUORS column
  
  all_cells$NAME <- as.numeric(all_cells$biodiverse_phylo_clusterS2)
}

all_cells$ELEMENT <- as.matrix(paste0(all_cells$x_coord, '_', all_cells$y_coord))
# View(all_cells)


clusters <- unique(all_cells$NAME, incomparables = FALSE)
clusters <- sort(clusters, decreasing = F)

groups <- list()
# print(groups)
for (i in 1:length(clusters)) {
  cluster <- all_cells[all_cells$NAME==clusters[i],]
  groups[[i]] <- cluster$ELEMENT
}
# View(groups[5])
# groups <- list(groups)
library(RColorBrewer)

library("ggplot2")
library("ggtree")

coloredtree <- groupOTU(tree, groups)
library("colorspace")

myPalette <- colorRampPalette(rev(brewer.pal(length(clusters), "Spectral")))
sc <- scale_colour_manual(values = c("black", myPalette(length(clusters))) )

myTree <- ggtree(coloredtree, aes(color=group)) + sc + theme(
  panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
#   panel.grid.minor = theme_blank(), 
#   panel.grid.major = theme_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA)
)

library("Cairo")
output_png <- paste0(taxon_group, '_', dendro_type, '_dendro.png')
output_pdf <- paste0(taxon_group, '_', dendro_type, '_dendro.pdf')
print('plotting at ')
print(paste0(output_folder,'/', output_png))

CairoPNG(width = 3000, height = 3000, file = paste(output_folder,'/', output_png, sep=""), bg = "transparent", units="px", dpi=300, title = "R Graphics Output") #
print(myTree)
dev.off()

CairoPDF(width = 7.5, height = 7.5, file = paste(output_folder,'/', output_pdf, sep=""), pointsize=10, bg = "white", paper = "special", pagecentre=TRUE) #
print(myTree)
dev.off()

closeAllConnections()




