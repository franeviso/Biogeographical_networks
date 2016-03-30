



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

PlotMe <- function(Axis_0, Axis_1, NAME, heatmap, for_rast, taxon_group, taxon_dir, shapefile, outline=F, outline_data=F, showNumbers=T, summary=F, summary_data=F, cont_region=F) {
  
  output_folder <- paste0(taxon_dir, 'figures/')
  dir.create(output_folder, showWarnings = FALSE)
  output_png <- paste0(taxon_group, '_', NAME, ".png")
  output_pdf <- paste0(taxon_group, '_', NAME, ".pdf")
  if (outline == T){
    output_png <- paste0(taxon_group, '_', NAME,'_',outline_data, ".png")
    output_pdf <- paste0(taxon_group, '_', NAME,'_',outline_data, ".pdf")
  } else if (cont_region == T) {
    output_png <- paste0(taxon_group, '_', NAME, "_colormap.png")
    output_pdf <- paste0(taxon_group, '_', NAME, "_colormap.pdf")
  }
  
  to_gen_rast <- for_rast
  
  if (outline == T) {
    all_cells_outline <- subset(to_gen_rast,  select=c(Axis_0, Axis_1, outline_data)) # r generate a rater based on the COLUORS column
    colnames(all_cells_outline) <- c("Axis_0", "Axis_1", "NAME")
    all_cells_outline$NAME <- as.numeric(all_cells_outline$NAME)
  }
  all_cells <- subset(to_gen_rast,  select=c(Axis_0, Axis_1, NAME)) # r generate a rater based on the COLUORS column
  colnames(all_cells) <- c("Axis_0", "Axis_1", "NAME")
  all_cells$NAME <- as.numeric(all_cells$NAME)
  
  map_data   <- readShapePoly(shapefile)
  map_extent <- extent(map_data)
  #----------------------------------------------------------------------
  
  # creating cluster numbering
  if (heatmap == F || outline == T) {
    myData <- all_cells
    if (outline == T){
      myData <- all_cells_outline
    }
    main <- as.data.frame(myData[,1:3])
    # View(main)
    
    main_clusters <- unique(main[,3], incomparables = FALSE)
    main_clusters <- sort(main_clusters, decreasing = F)
    
    central_points <- matrix(,length(main_clusters),3) 
    print(length(main_clusters))
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
  }
  
  myFont <- choose_font(c("HelvLight", "Arial", "sans"), quiet = TRUE) #load a font if available
  
#   df <- data.frame()
  p1 <- ggplot(data=all_cells)+
    geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="black", fill="transparent") +
    coord_fixed() 
  
  if (showNumbers == T) {
    numsize = 10
  } else {
    numsize = 4
  }
  if (heatmap == F) {
    
    # Print the summaries
    if (summary != F) {
      
      species_list <- summary_data[!grepl("\\d",summary_data$node), ]
      
      locations_list <- for_rast
      # locations_list <- summary_data[grepl("p_\\d",summary_data$node), ]
      # View(locations_list)
      regions <- arrange_(distinct_(select_(species_list, NAME)), NAME)
      if (cont_region == T) {
        new_all_cells <- all_cells
        central_points$region_CWE <- NA
      }
      output_summary <- paste0(taxon_group, '_', NAME, "_summary.txt")
      sink(paste0(output_folder,'/', output_summary))
      
      for (region in regions[[NAME]]) {
        in_region <- as.data.frame(species_list[species_list[, NAME] == toString(region), ])
        in_region_locations <- as.data.frame(locations_list[locations_list[, NAME] == toString(region), ])
#         if (region == 1) {
#           View(in_region_locations)
#         }
        region_WE <- sum(1/in_region$modularity_conv_num_links)
        region_CWE <- region_WE/nrow(in_region)
#         region_WE <- sum(in_region$network_conv_degree)
#         region_CWE <- 1/region_WE/nrow(in_region)
        cat( paste(region, ", number of species: ", nrow(in_region), ", number of locations: ", nrow(in_region_locations),
                   ", relative diversity: ", nrow(in_region)/nrow(in_region_locations), ", WE: ", region_WE, ", CWE: ", region_CWE ,
                   ", species: ", noquote(toString(in_region$node)), sep=" "))
        cat(". \n")
        if (cont_region == T) {
          central_points[region, ]$region_CWE <- region_CWE 
          new_all_cells[all_cells[, 'NAME'] == toString(region),]$NAME <- region_CWE
        }
      }
      sink()
      # View(new_all_cells)
    }
    
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
    
    myPalette <- colorRampPalette(rev(brewer.pal(length(main_clusters), "Spectral")))
    #  bias = 2
    # myPalette2 <- colorRampPalette(brewer.pal(100, "YlOrRd"), bias = 2)
    myPalette2 <- colorRampPalette(rev(brewer.pal(100, "Spectral")))
    if (cont_region == F) {
      sf <- scale_fill_manual(values = myPalette(length(main_clusters)) )
      sc <- scale_colour_manual(values = myPalette(length(main_clusters)), guide='none' )
      p1 <- p1 +  geom_tile(aes(Axis_0, Axis_1, fill=factor(NAME)), alpha=0.50) 
    } else {
      sf <- scale_fill_gradientn(colours = myPalette2(100), name = "")
      sc <- scale_colour_gradientn(colours = myPalette2(100), guide='none' )
      # limits=c(0, 0.5)
      p1 <- p1 +  geom_tile(data=new_all_cells, aes(Axis_0, Axis_1, fill=NAME), alpha=0.80) 
    }
    p1 <- p1  + sc+ sf+  geom_path(data=outline_pol.df, color="black", aes(long,lat,group=group) )
    
    # if (cont_region == F){
      theta <- seq(pi/8, 2*pi, length.out=16)
      xo <- diff(range(central_points$x))/200
      yo <- diff(range(central_points$y))/200
      for(i in theta) {
        p1 <- p1 + geom_text( data=central_points,
                              bquote(aes(x=x+.(cos(i)*xo),y=y+.(sin(i)*yo),label=number)), 
                              size=numsize, colour='black', family = myFont, fontface = 'bold' )
      }
      # View(central_points)
      if (cont_region == F){
        p1 <- p1 +  geom_text(aes(x=x, y=y, label= number, colour=factor(number)), data=central_points, size=numsize, family = myFont, fontface = 'bold')
      } else {
        p1 <- p1 +  geom_text(aes(x=x, y=y, label= number, colour=region_CWE), data=central_points, size=numsize, family = myFont, fontface = 'bold')
      }
      p1 <- p1+ theme( legend.position = "none",
                                title = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont)
                              )
    # }
    # p1 <- p1 + ggtitle(NAME)
    if (cont_region == T){
      p1 <- p1 +  theme( 
        #         legend.key.height = unit(1.1, "cm"), legend.margin   = unit(2, "cm"),
        legend.key.width  = unit(2.2, "cm"), #legend.position = c(.5, 0.93),
        # legend.key.height = unit(2.2, "cm"), #legend.margin   = unit(2, "cm"),
        #         legend.key.width  = unit(6.2, "cm"), legend.position = c(.5, 0.93),
        #         legend.direction  ='horizontal',
        legend.position = "top",
        title = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont)
        # legend.title      = element_text(colour = 'black', angle = 00, size=rel(1), face = 'plain', family = myFont)
        #         legend.text       = element_text(colour = 'black', angle = 0, size=rel(2), face = 'plain', family = myFont)
      )
    }

    
    
  } else {
#     scale_color_gradient2(..., low = muted("red"), mid = "white", high = muted("blue"), midpoint = 0, space = "rgb", na.value = "grey50", guide = "colourbar")
    myPalette <- colorRampPalette(rev(brewer.pal(100, "Spectral")))
    # myPalette <- colorRampPalette(brewer.pal(100, "YlOrRd"), bias = 2)
    # sc <- scale_fill_gradientn(colours = myPalette(100), name = NAME)
    sc <- scale_fill_gradientn(colours = myPalette(100), name = "")
    # sc <- scale_fill_gradientn(colours = myPalette(100), name = NAME, guide = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust=0.5, title.vjust=0.9, label.position="bottom", label.hjust = 0.5, label.vjust = 0.5, raster=FALSE))
    p1 <- p1  + geom_tile(aes(Axis_0, Axis_1, fill=NAME), alpha=0.8)
    # p1 <- p1 + ggtitle(NAME)
    
    
    p1 <- p1 + sc  +  theme( 
#         legend.key.height = unit(1.1, "cm"), legend.margin   = unit(2, "cm"),
        # legend.key.width  = unit(3.2, "cm"), #legend.position = c(.5, 0.93),
        legend.key.width  = unit(2.2, "cm"), #legend.position = c(.5, 0.93),
        # legend.key.height = unit(2.2, "cm"), #legend.margin   = unit(2, "cm"),
#         legend.key.width  = unit(6.2, "cm"), legend.position = c(.5, 0.93),
#         legend.direction  ='horizontal',
        legend.position = "top",
        title = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont)
        # legend.title      = element_text(colour = 'black', angle = 00, size=rel(1), face = 'plain', family = myFont)
#         legend.text       = element_text(colour = 'black', angle = 0, size=rel(2), face = 'plain', family = myFont)
      )
    
    
    if (outline == T) {
      
      layer_shape <- substr(basename(shapefile), 1, (nchar(basename(shapefile))-4))
      shape <- readOGR(shapefile, layer=layer_shape) #read the australia outline shape file
      proj <- projection(shape)
      
      spg <- all_cells_outline
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
      
      myPalette <- colorRampPalette(rev(brewer.pal(length(main_clusters), "Spectral")))
      sf <- scale_fill_manual(values = myPalette(length(main_clusters)) )
      sc <- scale_colour_manual(values = myPalette(length(main_clusters)), guide="none" )
      
      # p1 <- p1 +  geom_tile(aes(Axis_0, Axis_1, fill=factor(NAME)), alpha=0.40) 
      p1 <- p1 +  geom_path(data=outline_pol.df, color="black", aes(long,lat,group=group) )
      
      
      if (showNumbers == T) {
        theta <- seq(pi/8, 2*pi, length.out=16)
        xo <- diff(range(central_points$x))/200
        yo <- diff(range(central_points$y))/200
        for(i in theta) {
          p1 <- p1 + geom_text( data=central_points,
            bquote(aes(x=x+.(cos(i)*xo),y=y+.(sin(i)*yo),label=number)), 
            size=numsize, colour='black', family = myFont, fontface = 'bold' )
        }
        
        p1 <- p1 + sc + geom_text(aes(x=x, y=y, label= number, colour=factor(number)), data=central_points, size=numsize, family = myFont, fontface = 'bold' )
      }
      # p1 <- p1 + sc+ sf+ theme( legend.position = "none",
                                # title = element_text(colour = 'black', angle = 0, size=rel(1), face = 'plain', family = myFont)
      # )
      # p1 <- p1 + ggtitle(paste0(NAME, ' ', outline_data))
    }
    
  }
  # ------------------------------------ saving stuff ----------------------------
  
#   p1 + geom_polygon(data=map_data, aes(x=long, y=lat, group = group),colour="black", fill="transparent")
#   p1 +  coord_fixed()
  p1 <- p1 + theme( 
      text = element_text(family = myFont),
      strip.background = element_blank(),
      axis.line    = element_blank(),      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),      axis.ticks   = element_blank(),
      axis.title.x = element_blank(),      axis.title.y = element_blank(),
      panel.grid        = element_blank(),panel.background=element_blank(),
      plot.background   = element_rect(colour = "white", fill="white", size = 1), plot.margin=unit(c(0,0,-0.61,-0.61), "line")
    )
  
  
  print('plotting at ')
  print(paste0(output_folder,'/', output_png))
  if (heatmap == T || cont_region == T) {
    width = 2323
    height = 2000
  } else {
    width = 2323
    height = 2486
  }
#   CairoPNG(width = width, height = height, file = paste(output_folder,'/', output_png, sep=""), canvas="white", bg = "white", units="px", dpi=300, title = "R Graphics Output") #
#   print(p1)
#   dev.off()
  
  if (heatmap == T || cont_region == T) {
    width = 5
    height = 5
  } else {
    width = 5
    height = 5
  }
  CairoPDF(width = width, height = height, file = paste(output_folder,'/', output_pdf, sep=""), pointsize=10, bg = "white", paper = "special", pagecentre=T) #
  
  print(p1)
  dev.off()
  closeAllConnections()
}

# loading data

# inputs...


taxon_group <- "acaciaPhylo"

# shapefile <- "Z:/!my_working/data/shape/coastline_lamberts.shp"
shapefile <- "Z:/!my_working/data/shape/coastline_albers.shp"




# taxon_group <- "acaciaPhylo"
taxon_dir <- paste0("Z:/!my_working/data/", taxon_group, "/")
for_rast <- paste0(taxon_dir, taxon_group, "_almalgamated_data_geo_only.csv")
for_rast_sum <- paste0(taxon_dir, taxon_group, "_almalgamated_data.csv")
# for_rast <- paste0(taxon_dir, taxon_group, "_tablegroup.csv")
# shapefile <- "Z:/!my_working/data/shape/coastline_lamberts.shp"
# shapefile <- "Z:/!my_working/data/shape/coastline_albers.shp"

# shapefile <- choose.files (caption="select shape file", default="Z:/!my_working/data/shape/file")

for_rast <- read.csv(for_rast, header=T,sep=",")
for_rast_sum <- read.csv(for_rast_sum, header=T,sep=",")
# for_rast$network_CWE <- for_rast$map_equation_rw_flow / for_rast$modularity_num_links

# View(for_rast)
# calling the functions

#for paper maps...

PlotMe('x_coord', 'y_coord', 'modularity_conv_moduleID', F, for_rast, taxon_group, taxon_dir, shapefile, summary=T, summary_data = for_rast_sum)
# PlotMe('x_coord', 'y_coord', 'modularity_phylo_moduleID', F, for_rast, taxon_group, taxon_dir, shapefile, summary=T, summary_data = for_rast_sum)
# 
# 
PlotMe('x_coord', 'y_coord', 'biodiverse_CWE', T, for_rast, taxon_group, taxon_dir, shapefile)
# 
# PlotMe('x_coord', 'y_coord', 'modularity_conv_participation_coefficient', T, for_rast, taxon_group, taxon_dir, shapefile)
# PlotMe('x_coord', 'y_coord', 'modularity_phylo_participation_coefficient', T, for_rast, taxon_group, taxon_dir, shapefile)
# 
# 
PlotMe('x_coord', 'y_coord', 'map_equation_conv_moduleID', F, for_rast, taxon_group, taxon_dir, shapefile, summary=T, summary_data = for_rast_sum)
# PlotMe('x_coord', 'y_coord', 'map_equation_phylo_moduleID', F, for_rast, taxon_group, taxon_dir, shapefile, summary=T, summary_data = for_rast_sum)
# 
# PlotMe('x_coord', 'y_coord', 'map_equation_rw_moduleID', F, for_rast, taxon_group, taxon_dir, shapefile, summary=T, summary_data = for_rast_sum, showNumbers=F)
PlotMe('x_coord', 'y_coord', 'map_equation_rw_moduleID', F, for_rast, taxon_group, taxon_dir, shapefile, summary=T, summary_data = for_rast_sum, showNumbers=F, cont_region = T)
# PlotMe('x_coord', 'y_coord', 'map_equation_prw_moduleID', F, for_rast, taxon_group, taxon_dir, shapefile, summary=T, summary_data = for_rast_sum, showNumbers=F)
# 
PlotMe('x_coord', 'y_coord', 'biodiverse_clusterS2', F, for_rast, taxon_group, taxon_dir, shapefile)
# PlotMe('x_coord', 'y_coord', 'biodiverse_phylo_clusterS2', F, for_rast, taxon_group, taxon_dir, shapefile)
# 
# 
PlotMe('x_coord', 'y_coord', 'modularity_conv_participation_coefficient', T, for_rast, taxon_group, 
       taxon_dir, shapefile, outline = T, outline_data = 'modularity_conv_moduleID')
PlotMe('x_coord', 'y_coord', 'modularity_conv_participation_coefficient', T, for_rast, taxon_group, 
       taxon_dir, shapefile, outline = T, outline_data = 'biodiverse_clusterS2')
PlotMe('x_coord', 'y_coord', 'modularity_conv_participation_coefficient', T, for_rast, taxon_group, 
       taxon_dir, shapefile, outline = T, outline_data = 'map_equation_conv_moduleID')


# For the supporting information

PlotMe('x_coord', 'y_coord', 'biodiverse_CWE', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'biodiverse_WE', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'biodiverse_Richness', T, for_rast, taxon_group, taxon_dir, shapefile)

PlotMe('x_coord', 'y_coord', 'modularity_conv_participation_coefficient', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'modularity_conv_within_module_degree', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'map_equation_conv_flow', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'map_equation_rw_flow', T, for_rast, taxon_group, taxon_dir, shapefile)

PlotMe('x_coord', 'y_coord', 'network_conv_betweenness', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'network_conv_eigen_centrality', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'network_conv_alpha_centrality', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'network_conv_closeness', T, for_rast, taxon_group, taxon_dir, shapefile)

PlotMe('x_coord', 'y_coord', 'network_rw_betweenness', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'network_rw_eigen_centrality', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'network_rw_alpha_centrality', T, for_rast, taxon_group, taxon_dir, shapefile)
PlotMe('x_coord', 'y_coord', 'network_rw_closeness', T, for_rast, taxon_group, taxon_dir, shapefile)


