
library(ggplot2)
require(lattice)
require(RColorBrewer)
require(Cairo)
library(reshape2)




# functions

# plot one measure against another.
PlotMe <- function(data, plotme1, plotme2, output_folder, taxon_group) {
  
  plotting <- ggplot(data=data) + 
    #geom_point(data=random_data, aes(ENDW_RICHNESS,PE_WE_P, colour="RANDOMISED")) +
    geom_point(data=data, aes_string(plotme1, plotme2)) +
    stat_smooth(data=data, aes_string(plotme1, plotme2), method="lm", se=FALSE) +
    #labs(title = taxon)+
    #annotate("text", x = 4, y = -0.001, label = taxon) +
    theme(
      text = element_text(family = "HelvLight"),
      panel.background=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
    )+
    theme(
      panel.grid.major= element_line(colour = "grey"),
      panel.border = element_rect( fill="transparent", colour = "black"))
  correlation <- cor(data[[plotme1]], data[[plotme2]])
  plotting <- plotting + ggtitle(paste("correlation between", plotme1, "and", plotme2, "is ", correlation, sep=" "))
  
  output_png <- paste(taxon_group, plotme1, "vs", plotme2, ".png ", sep="_")
  print('plotting at ')
  print(paste0(output_folder,'/', output_png))
  CairoPNG(width = 2323, height = 2486, file = paste(output_folder,'/', output_png, sep=""), canvas="white", bg = "white", units="px", dpi=300, title = "R Graphics Output") #
  print(plotting)
  dev.off()
  
}

# metric for the distance between two sets.num_elem(A \cap B)/min_num_elem(A or B)
ClusterCompare <- function(myData, cluster1, cluster2, output_folder, taxon_group) {

#   myData <- data
#   cluster1 <- "biodiverse_clusterS2"
#   cluster2 <- "map_equation_moduleID"
  main <- as.data.frame(myData[[cluster1]])
  main$ID<-seq.int(nrow(main))
  # View(main)
  compare <- as.data.frame(myData[[cluster2]])
  compare$ID<-seq.int(nrow(compare))
  
  main_clusters <- unique(main[,1], incomparables = FALSE)
  main_clusters <- sort(main_clusters, decreasing = F)
  compare_clusters <- unique(compare[,1], incomparables = FALSE)
  compare_clusters <- sort(compare_clusters, decreasing = F)
  
  comparison_matrix <- matrix(,length(main_clusters),length(compare_clusters)) 
  for (i in 1:length(main_clusters)) {
    # print(i)
    main <- as.data.frame(main)
    module_members <- main[main[, 1] == main_clusters[i],]
    for (j in 1:length(compare_clusters)) {
      # print(j)
      compare_members <- compare[compare[, 1] == compare_clusters[j],]
      # value <- length(intersect(module_members[,2], compare_members[,2]))/min(length(module_members[,2]), length(compare_members[,2]))
      value <- length(intersect(module_members[,2], compare_members[,2]))/ length(compare_members[,2])
      comparison_matrix[i,j] <- value
      # print(value)
      
    }
    
  }
  
  rgb.palette <- colorRampPalette(c("white", "green"), space = "rgb")
  levelplot(t(comparison_matrix), col.regions=rgb.palette(120))
  # View(comparison_matrix)
  myPalette <- colorRampPalette(rev(brewer.pal(100, "YlGn")))
  sc <- scale_fill_gradientn(colours = rev(myPalette(100)), name = "Containment", limits = c(0, 1))
  
  df <- melt(comparison_matrix, varnames = c(cluster1, cluster2))
  
  
  df_box <- df
  df_box[[cluster1]] <- factor(df[[cluster1]])
  row_sub = apply(df_box[, -1], 1, function(row) !any(row <= 0.0001 ))
  ##Subset as usual
  df_box<-df_box[row_sub,]
  
  boxplotting <- ggplot(df_box, aes_string(cluster1, "value")) + geom_boxplot()+ scale_y_continuous(limits=c(0, 1))
  
  p1 <- ggplot(df, aes_string(cluster1, cluster2)) + geom_tile(aes(fill = value)) + sc
  p1 <- p1 + scale_y_continuous(breaks=seq(1,length(compare_clusters),1)) + scale_x_continuous(breaks=seq(1,length(main_clusters),1))
  output_png <- paste(taxon_group, cluster1, "vs", cluster2, ".png ", sep="_")
  print('plotting at ')
  print(paste0(output_folder,'/', output_png))
  # CairoPDF(width = 20, height = 60, file = paste(output_folder,'/', output_png, sep=""), pointsize=40, bg = "white", title = "R Graphics Output", version = "1.7", paper = "special", pagecentre=TRUE) #
  CairoPNG(width = 2000, height = 2000, file = paste(output_folder,'/', output_png, sep=""),  canvas="white", bg = "white", units="px", dpi=300, title = "R Graphics Output") #
  print(p1)
  dev.off()
  
  output_png <- paste(taxon_group, cluster1, "vs", cluster2, "boxplot.png ", sep="_")
  print('plotting at ')
  print(paste0(output_folder,'/', output_png))
  # CairoPDF(width = 20, height = 60, file = paste(output_folder,'/', output_png, sep=""), pointsize=40, bg = "white", title = "R Graphics Output", version = "1.7", paper = "special", pagecentre=TRUE) #
  CairoPNG(width = 2000, height = 2000, file = paste(output_folder,'/', output_png, sep=""),  canvas="white", bg = "white", units="px", dpi=300, title = "R Graphics Output") #
  print(boxplotting)
  dev.off()
}





# loading data


taxon_group <- "acacia100"

taxon_dir <- paste0("Z:/!my_working/data/", taxon_group, "/")

setwd(paste0(taxon_dir))
output_folder <- "analysis/"
dir.create(output_folder, showWarnings = FALSE)
data <- paste0(taxon_group, "_almalgamated_data_geo_only.csv")
# data <- paste0(taxon_group, "_almalgamated_data.csv")
data <- read.csv(data, header=T)
# plot(data)
data$network_CWE <- data$map_equation_rw_flow / data$modularity_num_links
# data$exp_CWE <- exp(data$biodiverse_CWE)
# data$maybe_CWE <- data$biodiverse_WE / data$modularity_num_links
# View(data)
# PlotMe(data, "exp_CWE", "network_CWE")
# # plotting data against eachother
# 
# PlotMe(data, "biodiverse_CWE", "maybe_CWE")
# 
# # 
# PlotMe(data, "biodiverse_CWE", "network_CWE", output_folder, taxon_group)
# # # PlotMe(data, "biodiverse_CWE", "modularity_within_module_degree")
# # # 
# PlotMe(data, "biodiverse_WE", "modularity_within_module_degree", output_folder, taxon_group)
# PlotMe(data, "biodiverse_WE", "map_equation_flow", output_folder, taxon_group)
# # # 
# PlotMe(data, "biodiverse_WE", "network_betweenness", output_folder, taxon_group)
# # PlotMe(data, "biodiverse_WE", "network_alpha_centrality", output_folder, taxon_group)
# # PlotMe(data, "biodiverse_WE", "network_eigen_centrality", output_folder, taxon_group)
# # PlotMe(data, "biodiverse_WE", "network_closeness", output_folder, taxon_group)
# # 
ClusterCompare(data, "modularity_moduleID", "map_equation_moduleID", output_folder, taxon_group)
ClusterCompare(data, "modularity_moduleID", "biodiverse_clusterS2", output_folder, taxon_group)

PlotMe(data, "modularity_participation_coefficient", "modularity_within_module_degree", output_folder, taxon_group)
# PlotMe(data, "modularity_participation_coefficient", "modularity_within_module_degree", output_folder, taxon_group)

# range weighted
PlotMe(data, "biodiverse_CWE", "network_CWE", output_folder, taxon_group)
PlotMe(data, "biodiverse_WE", "map_equation_rw_flow", output_folder, taxon_group)

PlotMe(data, "modularity_num_links", "map_equation_flow", output_folder, taxon_group)
PlotMe(data, "modularity_num_links", "network_betweenness", output_folder, taxon_group)
PlotMe(data, "modularity_num_links", "modularity_within_module_degree", output_folder, taxon_group)

# Phylo

PlotMe(data, "map_equation_phylo_flow", "biodiverse_PD", output_folder, taxon_group)
PlotMe(data, "network_phylo_betweenness", "biodiverse_PE_WE", output_folder, taxon_group)
PlotMe(data, "network_phylo_eigen_centrality", "biodiverse_PE_WE", output_folder, taxon_group)


# doing NMDS analysis
# vegan... NMDS analysis