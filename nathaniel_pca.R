library(vegan)
library(Cairo)


taxon_group <- "eucalypts_2016Phylo"
taxon_dir <- paste0("Z:/!my_working/data/", taxon_group, "/")
data_file <- paste0(taxon_dir, taxon_group, "_almalgamated_data_geo_only.csv")

data <- read.csv("/home/enc001/Documents/Nathaniel/data/eucalypts_2016Phylo/eucalypts_2016Phylo_almalgamated_data_geo_only.csv", header=T,sep=",") #read.csv(data_file, header=T,sep=",")



output_folder <- paste0(taxon_dir, 'pca/')
dir.create(output_folder, showWarnings = FALSE)

pca_variables <- data[,c("map_equation_conv_flow", 
                         "map_equation_rw_flow",
                         # "map_equation_phylo_flow", 
                         #"map_equation_prw_flow", 
                         "biodiverse_Richness", "biodiverse_WE", "biodiverse_CWE", 
                         # "biodiverse_PD", "biodiverse_PE_CWE", "biodiverse_PE_WE",
                         # "modularity_phylo_within_module_degree", "modularity_phylo_participation_coefficient", 
                         "modularity_conv_within_module_degree", "modularity_conv_participation_coefficient", 
                         "network_conv_alpha_centrality", "network_conv_eigen_centrality", "network_conv_betweenness", "network_conv_closeness",
                         "network_rw_alpha_centrality", "network_rw_eigen_centrality", "network_rw_betweenness", "network_rw_closeness"
                         # "network_phylo_alpha_centrality", "network_phylo_eigen_centrality", "network_phylo_betweenness", "network_phylo_closeness",
                         # "network_prw_alpha_centrality", "network_prw_eigen_centrality", "network_prw_betweenness", "network_prw_closeness"
                         )]

colnames(pca_variables) <- c("ME Flow", "ME RW Flow", "Richness", "WE", "CWE", "Modularity WMD", 
                             "Modularity PC", "Alpha-centrality", "Eigen-centrality", "Betweenness", "Closeness", 
                             "RW Alpha-centrality", "RW Eigen-centrality", "RW Betweenness", "RW Closeness")

stand.pca_variables <- scale(pca_variables)  

PCA <- rda(stand.pca_variables)

CairoPNG(width = 2000, height = 2000, file = paste(output_folder,'/conventional.png', sep=""), canvas="white", bg = "white", units="px", dpi=150, title = "R Graphics Output") #
biplot (PCA, display = 'species', type='points', scaling=3)
ordipointlabel(PCA, scaling = 3, display = c("species"), add=T, cex=1)
dev.off()

CairoPDF(width = 20, height = 20, file = paste(output_folder,'/conventional.pdf', sep=""), pointsize=35, bg = "white", paper = "special", pagecentre=TRUE) #
par(mar = c(4, 4, 1, 1))
biplot (PCA, display = 'species', type='points', scaling=3)
ordipointlabel(PCA, scaling = 3, display = c("species"), add=T, cex=1)
dev.off()

sink(paste(output_folder,'/conventional_PCA.txt', sep=""))
PCA
sink()

pca_variables_phylo <- data[,c(#"map_equation_conv_flow", 
                         "map_equation_phylo_flow",
                         "map_equation_prw_flow",
                         #"map_equation_prw_flow", 
                         # "biodiverse_Richness", "biodiverse_WE", "biodiverse_CWE", 
                         "biodiverse_PD", "biodiverse_PE_CWE", "biodiverse_PE_WE",
                         "modularity_phylo_within_module_degree", "modularity_phylo_participation_coefficient", 
                         # "modularity_conv_within_module_degree", "modularity_conv_participation_coefficient", 
                         # "network_conv_alpha_centrality", "network_conv_eigen_centrality", "network_conv_betweenness", "network_conv_closeness"
                         "network_phylo_alpha_centrality", "network_phylo_eigen_centrality", "network_phylo_betweenness", "network_phylo_closeness",
                         "network_prw_alpha_centrality", "network_prw_eigen_centrality", "network_prw_betweenness", "network_prw_closeness"
)]

stand.pca_variables_phylo <- scale(pca_variables_phylo)  

PCA_phylo <- rda(stand.pca_variables_phylo)

# biplot (PCA_phylo, display = 'species')
# plot(PCA_phylo, scaling = 3, type = "n")


CairoPNG(width = 2000, height = 2000, file = paste(output_folder,'/phylo.png', sep=""), canvas="white", bg = "white", units="px", dpi=150, title = "R Graphics Output") #
biplot (PCA_phylo, display = 'species', type='points', scaling=3)
ordipointlabel(PCA_phylo, scaling = 3, display = c("species"), add=T, cex=1)
dev.off()

CairoPDF(width = 15, height = 15, file = paste(output_folder,'/phylo.pdf', sep=""), pointsize=25, bg = "white", paper = "special", pagecentre=TRUE) #
biplot (PCA_phylo, display = 'species', type='points', scaling=3)
ordipointlabel(PCA_phylo, scaling = 3, display = c("species"), add=T, cex=1)
dev.off()


sink(paste(output_folder,'/phylo_PCA.txt', sep=""))
PCA
sink()

closeAllConnections()
