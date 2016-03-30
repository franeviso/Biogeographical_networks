


#function...
library(dplyr)
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


taxon_group <- "acaciaPhylo"

# shapefile <- "Z:/!my_working/data/shape/coastline_lamberts.shp"
shapefile <- "Z:/!my_working/data/shape/coastline_albers.shp"


# functions


outlineMe <- function(all_cells, shape) {
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
  return(outline_pol.df)
}

# inputs...


# blank theme
# new_theme_empty <- theme_bw()
# new_theme_empty$line <- element_blank()
# new_theme_empty$rect <- element_blank()
# new_theme_empty$strip.text <- element_blank()
# new_theme_empty$axis.text <- element_blank()
# new_theme_empty$plot.title <- element_blank()
# new_theme_empty$axis.title <- element_blank()
# new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")





taxon_dir <- paste0("Z:/!my_working/data/", taxon_group, "/")
for_rast <- paste0(taxon_dir, taxon_group, "_almalgamated_data_geo_only.csv")
for_rast_sum <- paste0(taxon_dir, taxon_group, "_almalgamated_data.csv")

for_rast <- read.csv(for_rast, header=T,sep=",")
for_rast_sum <- read.csv(for_rast_sum, header=T,sep=",")


output_folder <- paste0(taxon_dir, 'figures/')
dir.create(output_folder, showWarnings = FALSE)
output_png <- paste0(taxon_group, "_amalgamated.png")
output_pdf <- paste0(taxon_group, "_amalgamated.pdf")

to_gen_rast <- for_rast
Axis_0 = 'x_coord'
Axis_1 = 'y_coord'

map_data   <- readShapePoly(shapefile)
map_extent <- extent(map_data)
myFont <- choose_font(c("HelvLight", "Arial", "sans"), quiet = TRUE) #load a font if available

layer_shape <- substr(basename(shapefile), 1, (nchar(basename(shapefile))-4))
shape <- readOGR(shapefile, layer=layer_shape) #read the australia outline shape file

  
p1 <- ggplot()+
  geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="black", fill="transparent") +
  coord_fixed() 
  
  
for (count in 1:3) {
  if (count == 1) {
    NAME <- "map_equation_rw_moduleID"
  } else if (count == 2) {
    NAME <- "map_equation_conv_moduleID"
    
  } else if (count == 3) {
    NAME <- "modularity_conv_moduleID"
    
  }
  all_cells <- subset(to_gen_rast,  select=c(Axis_0, Axis_1, NAME)) # r generate a rater based on the COLUORS column
  colnames(all_cells) <- c("Axis_0", "Axis_1", "NAME")
  all_cells$NAME <- as.numeric(all_cells$NAME)
  
  #----------------------------------------------------------------------
  # numbering the points...
  myData <- all_cells
  main <- as.data.frame(myData[,1:3])
  # View(main)
  
  main_clusters <- unique(main[,3], incomparables = FALSE)
  main_clusters <- sort(main_clusters, decreasing = F)
  
  central_points <- matrix(,length(main_clusters),3) 
  for (i in 1:length(main_clusters)) {
    # print(i)
    main <- as.data.frame(main)
    module_members <- main[main[, 3] == main_clusters[i],]
    means <- colMeans(module_members)
    central_points[i, 1:2] <- means[1:2]
    central_points[i, 3] <- i
  }
  colnames(central_points) <- c("x",	"y",	"number")
  central_points <- as.data.frame(central_points)
  #----------------------------------------------------------------------
    
    if (count == 1) {
      numsize = 6
    } else if (count == 2) {
      numsize = 6
    } else if (count == 3) {
      numsize = 10
    }
    
    
    

    #---------------------------------------------------------------------------------
    if (count == 3) {
      myPalette <- colorRampPalette(rev(brewer.pal(length(main_clusters), "Spectral")))
      sc <- scale_colour_manual(values = myPalette(length(main_clusters)), guide="none" )
      outline_pol.df <- outlineMe(all_cells, shape)
      p1 <- p1 +  geom_path(data=outline_pol.df, color="black", aes(long,lat,group=group))
  
      
      theta <- seq(pi/8, 2*pi, length.out=16)
      xo <- diff(range(central_points$x))/200
      yo <- diff(range(central_points$y))/200
      for(i in theta) {
        p1 <- p1 + geom_text( data=central_points,
                              bquote(aes(x=x+.(cos(i)*xo),y=y+.(sin(i)*yo),label=number)), 
                              size=numsize, colour='black', family = myFont, fontface = 'bold' )
      }
      p1 <- p1 + sc+  geom_text(aes(x=x, y=y, label= number, colour=factor(number)), data=central_points, size=numsize, family = myFont, fontface = 'bold')
      
      # print(p1)
    #---------------------------------------------------------------------------------
    } else if (count == 1) {
      myPalette <- colorRampPalette(rev(brewer.pal(length(main_clusters), "Spectral")))
      sf <- scale_fill_manual(values = myPalette(length(main_clusters)), guide="none" )
      outline_pol.df <- outlineMe(all_cells, shape)
      p1 <- p1 +sf+  geom_tile(data= all_cells, aes(Axis_0, Axis_1, fill=factor(NAME)), alpha=0.5) 
      # p1 <- p1 +sf+   geom_path(data=outline_pol.df, color="grey", aes(long,lat,group=group))
      
#       
#       theta <- seq(pi/8, 2*pi, length.out=16)
#       xo <- diff(range(central_points$x))/200
#       yo <- diff(range(central_points$y))/200
#       for(i in theta) {
#         p1 <- p1 + geom_text( data=central_points,
#                               bquote(aes(x=x+.(cos(i)*xo),y=y+.(sin(i)*yo),label=number)), 
#                               size=numsize, colour='black', family = myFont, fontface = 'bold' )
#       }
#       p1 <- p1 +  geom_text(aes(x=x, y=y, label= number, fill=factor(number)), data=central_points, size=numsize, family = myFont, fontface = 'bold')
#       
      
      print(p1)
    
    #---------------------------------------------------------------------------------
    
    
    } else if (count == 2) {
      # all_cells <- as.data.frame(all_cells)
      sub_cells <- all_cells[ all_cells[,"NAME"] >= 6, ]
      outline_pol.df <- outlineMe(sub_cells, shape)
      # p1 <- p1 +  geom_tile(data= sub_cells, aes(Axis_0, Axis_1), color="red", alpha=0.1) 
      p1 <- p1 +  geom_path(data=outline_pol.df, color="red", aes(long,lat,group=group), size=1.2, alpha=1)
      
    }
    
    #---------------------------------------------------------------------------------
  
}

p1 <- p1  + theme( legend.position = "none",
                          title = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont)
)
# p1 <- p1 + ggtitle(NAME)
p1 <- p1 + theme( 
  text = element_text(family = myFont),
  strip.background = element_blank(),
  axis.line    = element_blank(),      axis.text.x  = element_blank(),
  axis.text.y  = element_blank(),      axis.ticks   = element_blank(),
  axis.title.x = element_blank(),      axis.title.y = element_blank(),
  panel.grid        = element_blank(),panel.background=element_blank(),
  plot.background   = element_rect(colour = "white", fill="white", size = 1), plot.margin=unit(c(0,0,-0.61,-0.61), "line")
) 

# ------------------------------------ saving stuff ----------------------------



print('plotting at ')
print(paste0(output_folder,'/', output_png))

width = 2323
height = 2486

# print(p1)

# grid.newpage()
# pushViewport( viewport( layout = grid.layout( 1 , 1 , widths = unit( 1 , "npc" ) ) ) ) 
# print( plot1, vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
# print( plot2, vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
# print( plot3, vp = viewport( layout.pos.row = 1 , layout.pos.col = 1 ) )
CairoPNG(width = width, height = height, file = paste(output_folder,'/', output_png, sep=""), canvas="white", bg = "white", units="px", dpi=300, title = "R Graphics Output") #
print(p1)
dev.off()
CairoPDF(width = 15, height = 15, file = paste(output_folder,'/', output_pdf, sep=""), pointsize=25, bg = "white", paper = "special", pagecentre=TRUE) #
mar = c(1, 1, 1, 1)
print(p1)
dev.off()
closeAllConnections()

