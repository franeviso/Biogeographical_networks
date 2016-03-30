
library(raster)
library(rgdal)
library(SDMTools)
library(extrafont)
library(grid)
library(gridExtra)
library(Cairo)
library(RColorBrewer)
library(ggplot2)
require("plyr")
require("rgdal") # requires sp, will use proj.4 if installed
require("maptools")

NAME <- 'biodiverse_phylo_clusterS2'
taxon_group <- "acaciaPhylo"
originalset <- 'Point_distribution_Australian_Phylogenetic_Diversity_Acacia.csv'
tree_file <- "acaciaPhylo_clean_tree_figtree.nwk"
coord1 <- "x_metres_EPSG_3577_Albers_Equal_Area"
coord2 <- "y_metres_EPSG_3577_Albers_Equal_Area"
taxa_label <- "Species"

taxon_dir <- paste0("Z:/!my_working/data/", taxon_group, "/")
tree_file_path <- paste0(taxon_dir, tree_file)
data_originalset <- paste0(taxon_dir, originalset)
datafile <- paste0(taxon_dir, taxon_group, "_almalgamated_data_geo_only.csv")
# for_rast <- paste0(taxon_dir, taxon_group, "_tablegroup.csv")
# shapefile <- "Z:/!my_working/data/shape/coastline_lamberts.shp"
shapefile <- "Z:/!my_working/data/shape/coastline_albers.shp"




for_rast <- read.csv(datafile, header=T,sep=",")


Axis_0 <- 'x_coord'
Axis_1 <- 'y_coord'
to_gen_rast <- for_rast
all_cells <- subset(to_gen_rast,  select=c(Axis_0, Axis_1, NAME))
colnames(all_cells) <- c("Axis_0", "Axis_1", "NAME")
all_cells$NAME <- as.numeric(all_cells$NAME)

map_data   <- readShapePoly(shapefile)

layer_shape <- substr(basename(shapefile), 1, (nchar(basename(shapefile))-4))
shape <- readOGR(shapefile, layer=layer_shape) #read the australia outline shape file
proj <- projection(shape)

spg <- all_cells
coordinates(spg) <- ~ Axis_0 + Axis_1
# coerce to SpatialPixelsDataFrame
gridded(spg) <- TRUE
# coerce to raster
rast <- raster(spg)
outline_raster <- rast
projection(outline_raster) <- CRS(proj)  #set its projection to the same as the shape file loaded above
outline_pol <- rasterToPolygons(outline_raster,n=4, dissolve=TRUE)
outline_pol@data$id = rownames(outline_pol@data)
outline_pol.points = fortify(outline_pol, region="id")
outline_pol.df = join(outline_pol.points, outline_pol@data, by="id")

ggplot()+
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="black", fill="transparent") +
  coord_fixed() + geom_path(data=outline_pol.df, color="black", aes(long,lat,group=group) )

# save polygons

dir.create(paste0(taxon_dir, "shpfiles"), showWarnings = FALSE)
setwd(paste0(taxon_dir, "shpfiles"))
writeOGR(obj=outline_pol, dsn=NAME, layer="NAME", driver="ESRI Shapefile")


# read shape file

plotpoly <- readShapePoly(paste0(NAME,"/NAME.shp"))
ggplot()+
  geom_polygon(data=plotpoly, aes(x=long, y=lat, group = group),colour="black", fill="transparent")+
  coord_fixed()

#-------------------------------------------
# loading dataset and assigning regions
#--------------------------------------------

#this file loads a spatial points dataset (csv file) and attaches values to the points from a shape file
library(sp)
library(rgdal)
library(maps)
library(raster)

file_loaded <- read.csv(file=data_originalset,header=T,sep=",")# load file
# View(file_loaded)
file_sub <- file_loaded
#View(file_sub)
coordinates(file_sub) <- c(coord1, coord2)# convert it to coordinates

# plot(file_sub[1,])
#wgs84proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #WGS84 proj4 projection
shape <- readOGR(paste0(NAME,"/NAME.shp"), "NAME") #read in a shape file to

projection(shape)
proj4string(file_sub) <- projection(shape)
#proj4string(file_sub) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # set the projections 
inside.shape <- !is.na(over(file_sub, as(shape, "SpatialPolygons"))) # get the points that are inside the polygon

mean(inside.shape) # what fraction of points are inside the polygon
#summary(shape)

file_sub$region <- over(file_sub, shape)$NAME # add a column to the dataset and populate it with a value from the polygons in the shape file
#file_clipped <- file_sub[!is.na(file_sub$region),] # clip all the ones that have no value for the new column
#View(over(file_sub, shape))
#View(file_sub)

write.csv(file_sub, file=paste0(NAME,"/", NAME, "_withregions.csv"), row.names = FALSE)#write t

#-------------------------------------------
# creating biogeobears input files
#--------------------------------------------
getwd()
presence_list <- read.csv(file=paste0(NAME,"/", NAME, "_withregions.csv"),header=T,sep=",")# 
# presence_list <- as.data.frame(file_sub)
# file_sub[, taxa_label] <- as.data.frame(file_sub[, taxa_label])

unique_species <- unique(presence_list[[taxa_label]])
all_regions <- sort(unique(presence_list[["region"]]))
max_regions <- length(unique(presence_list[["region"]]))
base_string <- ''
for (i in 1:max_regions) {
  base_string <-paste0('0', base_string)
}
# substring(base_string, 2, 2) <- '1'

occurance_file <- matrix(0, ncol = 2, nrow = length(unique_species))
# unique_species
i= 1
for (specie in unique_species) {
  species_subset <- presence_list[presence_list[, taxa_label] == specie,]
  unique_regions <- unique(species_subset[["region"]])
  occurance_file[i, 1] <- specie
  region_string <- base_string
  for (number in unique_regions) {
    region_num <- which(all_regions == number)
    substring(region_string, region_num, region_num) <- '1'
  }
  occurance_file[i, 2] <- region_string
  i = i+1
}

#-------------------------------------------
# trimming
#--------------------------------------------

# Now need to trim tree and datafile against each other.... Or shouldn't need to in this case?

library(tidyr)
library(ape)
library(phytools)
# setwd("C:/GIS-Datasets/marsupials_final")

tree <- read.newick(file=tree_file_path)
tree$tip.label
#tree <- read.tree(file=tree_file)

occurance_matrix_df <- as.data.frame(occurance_file)

trimmed_occurance_matrix <- occurance_matrix_df[occurance_matrix_df[["V1"]] %in% tree$tip.label, ]

# View(trimmed_occurance_matrix)
# deleted_taxa_occurance_matrix <- occurance_matrix_df[ ! occurance_matrix_df$Taxa %in% tree$tip.label, ]
# View(deleted_taxa_occurance_matrix)


number_taxa <- nrow(trimmed_occurance_matrix)
number_regions <- max_regions

phylip_format <- trimmed_occurance_matrix
phylip_format <- unite(phylip_format, merged, 2:ncol(trimmed_occurance_matrix), sep="")
#View(phylip_format)

regions <- noquote(seq.int(1, number_regions, 1))
  
#   
#   noquote(colnames(trimmed_occurance_matrix)) #get regions as list for first line of phylip (biogeobears needs this)
# regions <- noquote(regions[-1]) #remove the first value as it is not a region

# setwd("C:/GIS-Datasets/marsupials_final/biogeobears_supercomputer")

sink(paste0(NAME,"/", taxon_group, '_', NAME, ".phy"))
cat( paste(number_taxa, number_regions, sep=" "))
cat(" (")
cat(paste0(regions))
cat(") \n")
sink()
write.table(phylip_format, file=paste0(NAME,"/", taxon_group, '_', NAME, ".phy"),  row.names = FALSE, col.names=FALSE, append=TRUE, quote=FALSE, sep="\t")#write phylip file to input into biogeobears
write.tree(tree, file = paste0(NAME,"/", tree_file), append = FALSE, digits = 10, tree.names = FALSE)

