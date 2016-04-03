# Biogeographical_networks
Set of R scripts to construct and analyse biogeographical networks to detect bioregions


## How to Use these Scripts

This is a breakdown of how I did all of the analyses. All of the scripts are name spaced, so if you keep to the same naming structure things will be a breeze and you will only have to change the taxon_group at the top of each script.

### Running the Core Analysis

1. Load data into biodiverse. Load the tree into biodiverse. Get the label names to match if required (with the label remap in biodiverse), and then trim tree to basedata and basedata to tree, if using a tree. Run clustering and spatial analysis.
2. Export the basedata as a table grouped csv file, called "foldername\_label\_export.csv". Export the tree as a newick file. Open it up in Figtree and then save it as "foldername\_clean\_tree\_figtree.nwk" - this deals with biodiverse's slightly buggy tree output.

    *If conducting phylogenetic analysis (requires tree):*
3. **nathaniel\_load\_tree\_and\_biodiverse\_label\_export\_trim\_then\_prepare\_for\_generation\_of\_network\_files.R**: Prepares the files for the creation of the phylogenetic and range-weighted networks.
3. **nathaniel\_load\_files\_and\_create\_network\_matrix.R**: Creates the phylogenetic and phylogenetic range-weighted networks.

    *If not conducting phylogenetic analysis (no tree required):*
3. **nathaniel\_load\_files\_and\_create\_tips\_only\_network\_matrix.R**: Creates the tips-only unweighted and tips-only range-weighted networks.

    *Running the analysis:*
4. **nathaniel\_open\_cmd\_run\_analyses.R**: Calls the command line for either modularity or the map equation. Can run all network analyses from the command line using this script.

### Helper scripts

These are a set of scripts I made to make dealing with the results a little easier. They all require the following script to be run first

+ **nathaniel\_process\_results.R**: Gathers all of the available data and puts it into two csv files, one which contains only the geographic nodes (\_amalgamated\_data\_geo\_only), and the other which contains all nodes (\_amalgamated\_data). Need to select which data to collate at the start of the file. It also can add the biodiverse results - see the script itself for the naming conventions.

Once the above script is run, these other scripts acomplish everything else

#### Plotting
+ **nathaniel\_plot\_cluster\_map.r**: Plots the data onto a map over Australia. Can pick which variable, how it should be plotted and give summaries of the bioregions.
+ **nathaniel\_plot\_cluster\_map\_overlay.r**: Plots different bioregionalisations onto the same figure. Combines the range-weighted and unweighted map equation, and the modularity bioregionalisations.
+ **nathaniel\_plot\_cluster\_dendogram.r**: Plots the dendrograms from Biodiverse. Need to export the dendrogram as a newick file from Biodiverse first.
+ **nathaniel\_analyse\_results.r**: Outputs correlation plots for two variables, and also has graphing methods (box-plot and heatmap matrix) for comparing the bioregionalisation between two different methods.
+ **nathaniel\_pca.r**: Creates a PCA plot from chosen variables.

#### Miscellaneous
+ **nathaniel\_plot\_cluster\_map\_table.r**: Creates the latex tables for the SI, which summarise the number of modules in each bioregionalisation, average inverse range of their species, number of locations and taxa, and also lists a few of their most indicative taxa.
+ **nathaniel\_species\_metrics.r**: Outline script to examine species which rank highly on network variables, such as participation coefficient or betweenness.

#### Biogeobears Helpers
+ **nathaniel\_create\_biogeobears.r**: Creates the regions file for input into biogeobears.

#### Gephi Helpers
+ **Gython Script.py**: Used in the Jython command line in Gephi (Scripting Plugin). The line needed to run it is commented in the top of the script. To use it, load the network and import the almalgamated data. You will need to change either ID.x or ID.y to ID, depending on if you are using the conventional or phylogenetic network. It may also be nessecary to delete the last row imported by gephi (as this is null and screws up the next bit in the conventional case). Then set the type of regionalisation you wish to color by in the script (cluster\_type), and you're ready to run it. If you need more colors, run **nathaniel\_generate\_colors\_for\_gython.r** and add the string and number of colors neeeded to the dictionary in the Gython script. Afterwards, run Force Atlas 3D, and if you need to flip the network use the geometricTransformation plugin with a homeothetic transform of -1 (its hard to find, but there is a link to it in the Gephi forum).
 + **nathaniel\_generate\_colors\_for\_gython.r**: See above.
