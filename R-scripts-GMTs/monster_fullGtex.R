# Part I: makes an expression matrix based on gene library and 3000 randomly chosen samples

# R script to download selected samples
# Copy code and run on a local machine to initiate download

# Check for dependencies and install if missing
packages <- c("rhdf5","downloader")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  print("Install required packages")
  source("https://bioconductor.org/biocLite.R")
  biocLite("rhdf5")
  install.packages("downloader", dependencies=T)
}
library("rhdf5")
library("downloader")

destination_file = "clean_gtex.rda"
extracted_expression_file = "sample_expression_matrix.tsv"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
  print("Downloading compressed gene expression matrix.")
  url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
  download(url, destination_file, mode="wb")
} else{
  print("Local file already exists.")
}

# the GTEX version of expression creation
expression = gtex

# we want to do quantile normalization later on
# install package
# source("https://bioconductor.org/biocLite.R")
# biocLite("preprocessCore")

# load package
# library(preprocessCore)

# get the index of the rowname of <NA> and remove it
toRm = which(is.na(rownames(expression)))
expression = expression[-toRm, ]

# do not accept the rows that are all zero
expression = expression[rowSums(!as.matrix(expression)) < ncol(expression), ]

# aggregate the duplicate rows by median and give the final expression matrix
agg = aggregate(expression, list(rownames(expression)), median)
# agg = as.data.table(expression, keep.rownames = T)
# agg = ddply(agg, .(rn), numcolwise(median))
rownames(agg) = agg[,1]

expression = as.matrix(agg[, -1])

# Print file
# write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE)
# print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))

cat("Dimensions of Expression: ", dim(expression), "\n")

# we want to do quantile normalization on the entire expression matrix
nEx = normalize.quantiles(expression)
rownames(nEx) = rownames(expression)
colnames(nEx) = colnames(expression)


# Part II: Make the expression matrix only the genes also in the gmt file. To do so, read the gmt file.

gmtNames = c("ChEA_2016.txt", "ENCODE_TF_2015.txt", "GO_Biological_Process_2017.txt", 
		"GO_Cellular_Component_2017.txt", "GO_Molecular_Function_2017.txt", "KEA_2015.txt", 
		"KEGG_2016.txt", "MGI_Mammalian_Phenotype_Level_4.txt")

# This script just extracts all the genes from the H_P_O file and removes duplicates. It puts them into a vector.
# it first reads the gmt file

gtexCoAvg = list()
gtexAUC = list()

# HERE BEGINS THE FOR LOOP THAT PRODUCES AUC AND average correlation matrices FOR EACH OF THE GMT FILES 
for (txtFile in gmtNames) {

  print(txtFile)

  gmt = readLines(txtFile)

  go = list()
  for(line in gmt){
    # splits line into the tabbed words
    sp = unlist(strsplit(line, "\t"))
    
    # term is the name of the current phenotype
    term = sp[1]
    
    # t ends up being all the current phenotype's genes, but it's currently empty
    t = c()
    
    # starts at 3 because that's where the genes begin (the name of the phenotype, an empty space, then first gene)
    for(m in 3:length(sp)){
      sp1 = unlist(strsplit( sp[m], ","))
      t = c(t, sp1[1])
    }
    # holds all the phenotypes and their genes in the list format, where the name of the phenotype is the category name
    # Each phenotype is a vector of its genes
    # taking the intersection of the names in expression and the names in this particular list element/phenotype
    go[[length(go)+1]] = intersect(rownames(expression), t) 
    names(go)[length(go)] = term
  }

  # because of the "intersect()", we now have only the unique genes in both expression library and H_P_O file
  gene_vector = unlist(go, use.names = FALSE)
  gene_vector = unique(gene_vector, incomparables = FALSE, fromLAST = FALSE, nmax = NA)

  # now let's normalize matrix of the gene library genes
  # We can do this outside of the loop - should be the same for all steps
  # norm_exprTxt = normalize.quantiles(expression[gene_vector,])
  # rownames(norm_exprTxt) = rownames(expression[gene_vector,])
  # norm_transposeTxt = t(norm_exprTxt)
  libEx = nEx[gene_vector, ]


  # BEGINNING OF AUC CODE
  # ok - let's loop now.
  steps = 20
  total = nrow(expression)
  size = total %/% steps

  # scales to a range
  scale_vector <- function(x, start, end)
    (x - min(x)) / max(x - min(x)) * (end - start) + start

  # a matrix for new GMT making and a vector for violin plots
  allCrrAvg = matrix(NA, nrow = nrow(expression), ncol = length(go))
  rownames(allCrrAvg) = rownames(expression)
  colnames(allCrrAvg) = names(go)
  auc = c()

  for (i in 1:steps) {
    
    # so instead of normalizing subsets (again), take only the subset you need...
    exprStep = nEx[(i + size * (i - 1) ):(min(i + size * i, total)),]
    
    # making the correlation between all genes and genes in the text file
    stepCorrelation = cor(x = t(exprStep), y = t(libEx))
    rownames(stepCorrelation) = rownames(exprStep)
    colnames(stepCorrelation) = rownames(libEx)

    stepCorrelation[which(stepCorrelation == 1)] = NaN
    # removes correlation between protein/gene and itself
    # for (j in 1:length(stepCorrelation)) {

    #   v1 = as.vector(arrayInd(j, dim(stepCorrelation)))

    #   if (rownames(correlation)[v1[1]] == colnames(stepCorrelation)[v1[2]]) {
    #     stepCorrelation[j] = NA
    #   }

    # }
    
    # now calculate the means of each gene's correlation to each of the gene sets
    co_means = matrix(NA, nrow = nrow(stepCorrelation), ncol = length(names(go)))
    rownames(co_means) = rownames(stepCorrelation)
    colnames(co_means) = names(go)
    
    # more efficient matrix making - taking a subset of genes based on gene set and averaging the correlations all at once
    for (k in 1:length(go)) {
      gene_set = go[[k]]
      
      if (length(gene_set) == 1) {
        co_means[, k] = stepCorrelation[, gene_set]
      } else {
        # now we subset and take the row means
        # co_sums[, j] <- rowSums(stepCorrelation[,gene_set], na.rm = TRUE)
        co_means[, k] = rowMeans(stepCorrelation[,gene_set], na.rm = TRUE)
      }
    }
    
    allCrrAvg[(i + size * (i - 1) ):(min(i + size * i, total)) , ] = co_means
    
    # so we construct the auc vector here - remember, the auc vector must be initialized outside of the giant step loop
    # now apply to all rows/genes
    for (current_gene in rownames(co_means)) {
      
      # this will eventually identify which sets have the gene and which sets do not
      containing_gene = co_means[current_gene,]
      
      # sorts the correlation averages
      containing_gene = sort(containing_gene, decreasing = TRUE)
      
      # for each gene_set in which the gene is a member, we want a 1 to be associated; for the rest, a zero.
      for (n in 1:length(go)) {
        if (current_gene %in% go[[n]]) {
          containing_gene[names(go)[n]] = 1
        } else {
          containing_gene[names(go)[n]] = 0
        }
      }
      
      cumulative = cumsum(containing_gene)
      
      # Part II: Now going to make a vector of AUC for all genes. 
      
      scaled_y = scale_vector(cumulative, 0, 1)
      scaled_x = scale_vector(1:length(cumulative), 0, 1)
      
      # find the area under the curve
      auc[current_gene] = trapz(scaled_x, scaled_y)
    }
    
  }

  # gtexCoAvg[[length(gtexCoAvg) + 1]] = allCrrAvg
  # names(gtexCoAvg)[length(gtexCoAvg)] = txtFile
  
  # gtexAUC[[length(gtexAUC) + 1]] = auc
  # names(gtexAUC)[length(gtexAUC)] = txtFile
	
  # create text files and do not store "allCrrAvg"
  # so all the new GMTs will be saved as text files
  toGmt_size(allCrrAvg, paste0("newGmtSize_GTEX_", txtFile), go)
  toGmt_z(allCrrAvg, paste0("newGmt(3)_GTEX_", txtFile), 3)

  write.csv(auc, paste0("tempGtexAUC", txtFile, ".csv"))

  print(mean(auc, na.rm = T))

}
