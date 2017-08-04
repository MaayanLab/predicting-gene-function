# test-site-files

## GMTs
The main idea is to predict gene function by creating new gene libraries based on their expression data correlations.

### Method
#### Creating the New GMTs
1. Download the expression matrix data set of interest and move it to your working directory.
    1. The ARCHS4 Human and Mouse data sets are available on the ARCHS4 website/chrome extension.
    2. The GTEX and CCLE data sets are available on their respective websites.
2. Run monster_fullQN.R for the Human and Mouse data sets. Run monster_fullGtex.R for GTEX, and run monster_fullCcle.R for CCLE.
3. This will create new GMT text files in your working directory.

#### Making Violin Plots
1. Follow the instructions within the commented R script, graphics(2).R. Use it as a template for creating new violin plots
    1. Note: You will need the AUC data for each library to make its violin plot. This data is available in the data folder.

## PPIs
The main idea is to predict new protein-protein interactions based on their expression data correlations.

### Method
#### Creating the New PPIs
1. Run/Source ppiConversionFunctions.R. 
2. Run/Source monster_fullPpi_steps(2).R. This will create the following:
    1. Text files of the new PPIs, formatted as sets.
    2. A .csv file of all areas under a ROC curve (for each gene in the original PPI).
    3. An R list with a filtered version of the original PPI, which is necessary for the creation of Venn diagrams.

#### Making Venn Diagrams
1. Convert the newly generated PPIs from set form to paired form (a function is available for this in the ppiConversionFunctions.R script.)
2. Use the pertinent Venn diagram R script.
    1. 

## Miscellaneous
1. In creating the violin plots, the violin data from the data folder should be uploaded beforehand.
2. For convenience, an R vector-ready list of the Enrichr library names is in the data folder. 
