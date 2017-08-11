# predicting-gene-function

## GMTs
The main idea is to predict gene function by creating new gene libraries based on their expression data correlations.

### Method
#### Creating the New GMTs
1. Download the expression matrix data set of interest and move it to your working directory.
    1. The ARCHS4 Human and Mouse data sets are available on the ARCHS4 website/chrome extension.
    2. The GTEX and CCLE data sets are available on their respective websites.
2. Download the library or libraries of interest and move them to your working directory (suggested: download from Enrichr. These libraries must be in the .gmt file format.).
3. Open the R console. Install and load Bioconductor package preprocessCore. Install and load R package pracma.
4. Ensure that functions have been created (specifically, the toGmt functions available in the R-Scripts-GMT folder). Run the toGmt scripts (only toGmt_size.R and toGmt_z.R are needed). 
5. Run monster_fullQN.R for the Human and Mouse data sets. Run monster_fullGtex.R for GTEX, and run monster_fullCcle.R for CCLE. Depending on how many libraries you are running and the size of the libraries, this can take anywhere from thirty minutes to one day.
    1. monster_fullQN.R: You will have to adjust the code, depending on whether you want to analyze the Human or Mouse data sets.
6. This will create new GMT text files in your working directory. This will also record all AUCs for each gene in the pertinent data set in relation to the pertinent library.

#### Making Violin Plots
1. Install the R packages ggplot2 and gridExtra.
2. Follow the instructions within the commented R script, graphics(2).R. Use it as a template for creating new violin plots
    1. Note: You will need the AUC data for each library to make its violin plot. This data is available in the data folder.
    2. This data needs to be converted to an R data.frame, and all NA values should be removed. The graphics(2).R script should provide indications as to how the violin data must be formatted before making violin plots with ggplot2.

##### Rough Reminders and Notes
* The graphics(2).R script is designed to take in a data.frame with columns "Genes", "AUC", and "Expression". The rows will be each gene name. This single data.frame holds all genes from each co-expression data set and makes a violin plot comparing their AUCs with respect to a single library.
* Remove the NA values. 
* Remember to convert from the uploaded .csv to an R data.frame. 
* Once you run the code for getting the graph, you can remove the legend at will by using get_legend (in the script) and then setting the legend.position to "none". For example:
```R
legend = get_legend(pCh)
pCh = pCh + theme(legend.position = "none")
```

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
1. Install and load the R packages data.table and VennDiagram.
2. Convert the newly generated PPIs from set form to paired form (a function is available for this in the ppiConversionFunctions.R script).
3. Use the pertinent Venn diagram R script.
    1. For example, a triple Venn diagram of all the old PPI files can be created by the R script, intVenn_allOld.R.

## Miscellaneous
1. In creating the violin plots, the violin data from the data folder should be uploaded beforehand.
2. For convenience, an R-vector-ready list of the Enrichr library names is in the data folder. 
3. In the R-scripts-GMTs folder, the R script tTest.R records t-test results about 8 Enrichr libraries and 4 expression data sets. 
4. In the gmtSite folder, the R scripts of idg_functions.R and idg_predictions.R will help make organized lists of predicted functions for each protein in the NIH's program, Illuminating the Druggable Genome. 
5. Overall, the packages used include:
    1. pracma (finding area under ROC curve)
    2. ggplot2 (plotting violin plots)
    3. gridExtra (plotting multiple violin plots at once)
    4. RColorBrewer (color palettes for violin plots)
    5. VennDiagram (making PPI Venn diagrams)
    6. preprocessCore (from Bioconductor; used in quantile normalization of expression data)
   
