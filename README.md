# predicting-gene-function

## GMTs
The main idea is to predict gene function by creating new gene libraries based on their expression data correlations.

### Method
#### Creating the New GMTs
1. Download the expression matrix data set of interest and move it to your working directory.
    1. The ARCHS4 Human and Mouse data sets are available on the ARCHS4 website/chrome extension.
    2. The GTEX and CCLE data sets are available on their respective websites.
2. Download Bioconductor package preprocessCore. Download R packages pracma and ggplot2 (for later violin plot creation).
3. Ensure that functions have been created (specifically, the toGmt functions available in the R-Scripts-GMT folder).
3. Run monster_fullQN.R for the Human and Mouse data sets. Run monster_fullGtex.R for GTEX, and run monster_fullCcle.R for CCLE.
4. This will create new GMT text files in your working directory.

#### Making Violin Plots
1. Install the R packages ggplot2 and gridExtra.
2. Follow the instructions within the commented R script, graphics(2).R. Use it as a template for creating new violin plots
    1. Note: You will need the AUC data for each library to make its violin plot. This data is available in the data folder.
    2. This data needs to be converted to an R data.frame. The graphics(2).R script should provide indications as to how the violin data must be formatted before making violin plots with ggplot2.

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
1. Install the R packages data.table and VennDiagram.
2. Convert the newly generated PPIs from set form to paired form (a function is available for this in the ppiConversionFunctions.R script).
3. Use the pertinent Venn diagram R script.
    1. For example, a triple Venn diagram of all the old PPI files can be created by the R script, intVenn_allOld.R.

## Miscellaneous
1. In creating the violin plots, the violin data from the data folder should be uploaded beforehand.
2. For convenience, an R-vector-ready list of the Enrichr library names is in the data folder. 
3. In the R-scripts-GMTs folder, the R script tTest.R records t-test results about 8 Enrichr libraries and 4 expression data sets. 
4. Overall, the packages used include:
    1. pracma (finding area under ROC curve)
    2. ggplot2 (plotting violin plots)
    3. gridExtra (plotting multiple violin plots at once)
    4. RColorBrewer (color palettes for violin plots)
    5. VennDiagram (making PPI Venn diagrams)
    6. preprocessCore (from Bioconductor; used in quantile normalization of expression data)
   
