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

destination_file = "human_matrix_download.h5"
# destination_file = "mouse_matrix.h5"
# destination_file = "clean_gtex.rda"
extracted_expression_file = "sample_expression_matrix.tsv"

# Check if gene expression file was already downloaded, if not in current directory download file 
# form repository
if(!file.exists(destination_file)){
  print("Downloading compressed gene expression matrix.")
  url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
  download(url, destination_file, mode="wb")
  } else{
    print("Local file already exists.")
  }

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/samples")
tissue = h5read(destination_file, "meta/tissue")
# library = h5read(destination_file, "meta/library")
genes = h5read(destination_file, "meta/genes")
# instrument = h5read(destination_file, "meta/instrument")

# Identify columns to be extracted
sample_locations = sample(1:length(samples), 3000)

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
# rownames(expression) = genes
# colnames(expression) = samples[sample_locations]

# the GTEX version of expression creation
# expression = gtex
# rownames(expression) = genes
# colnames(expression) = samples

# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))

rownames(expression) = genes
colnames(expression) = samples[sample_locations]

cat("Dimensions of Expression:", dim(expression), "\n")

# we want to do quantile normalization on the entire expression matrix
nEx = normalize.quantiles(expression)
rownames(nEx) = genes
colnames(nEx) = samples[sample_locations]

# Part II: Make the expression matrix only the genes also in the gmt file. To do so, read the gmt file.

gmtNames = c("biogrid_ppi_2017_06_09.sig", "HuMAp_ppi_2017_06_07.sig", "bioplex_ppi_2017_06_07.sig")
fileNames = c("Biogrid", "Humap", "Bioplex")

# extracts all the genes from the H_P_O file and removes duplicates. It puts them into a vector.
# it first reads the gmt file

ppiOld = list()

scale_vector <- function(x, start, end)
(x - min(x)) / max(x - min(x)) * (end - start) + start

# HERE BEGINS THE FOR LOOP THAT PRODUCES AUC AND AVERAGE CORRELATION FOR EACH OF THE GMT FILES 
for (a in 1:length(gmtNames)) {

  print(gmtNames[a])
  gmt = readLines(gmtNames[a])

  v1 = c()
  v2 = c()

  for (c in 1:length(gmt)) {
    sp = unlist(strsplit(gmt[c], "\t"))

    v1[c] = sp[1]
    v2[c] = sp[6]

  }

  proteins = c(v1, v2)

  proteins = unique(proteins, incomparables = FALSE, fromLAST = FALSE, nmax = NA)
  proteins = intersect(rownames(expression), proteins)
  cat("length(proteins) =  ", length(proteins), "\n")

  # now let's normalize sub-matrix of expression, where we just take the "genes" that are the proteins 
  # in the ppi
  # norm_exprTxt = normalize.quantiles(expression[proteins,])
  # rownames(norm_exprTxt) = rownames(expression[proteins,])
  # norm_transposeTxt = t(norm_exprTxt)
  libEx = nEx[proteins, ]

  # read in the ppi and make sets
  go = vector("list", length(proteins))
  names(go) = proteins

  for (b in 1:length(gmt)) {

    sp = unlist(strsplit(gmt[b], "\t"))
    
    if (sp[1] != sp[6] && sp[1] %in% proteins && sp[6] %in% proteins) {
    	term = sp[1]
    	go[[term]] = c(go[[term]], sp[6])

    	term2 = sp[6]
    	go[[term2]] = c(go[[term2]], sp[1])
    }
  }

  # remove the null values
  go = go[vapply(go, Negate(is.null), NA)]
  # remove duplicates
  go = lapply(go, unique) 
  # remove the proteins with less than 5 interactions
  go = go[which(lengths(go) >= 5)]
  # intersect with expression
  # go = lapply(go, function(x) intersect(rownames(expression), x))
  cat("length(go) =  ", length(go), "\n")

  # BEGINNING OF AUC CODE

  # # removes correlation between protein/gene and itself
  # for (i in 1:length(correlation)) {
  #   v1 = as.vector(arrayInd(i, dim(correlation)))
  #   if (rownames(correlation)[v1[1]] == colnames(correlation)[v1[2]]) {
  #     correlation[i] = NA
  #   }
  # }

  # indices = which(correlation %in% sort(correlation, decreasing = T)[1:length(l1)])
  # for (j in 1:length(indices)) {
  #   pair = arrayInd(indices[j], dim(correlation))
  #   cat(rownames(correlation)[pair[1]], colnames(correlation)[pair[2]], file = fileNames[a], 
        # sep = "\t\t\t\t\t")
  # }

  # ok - let's loop now.
  steps = 20
  total = nrow(expression)
  size = total %/% steps

  # a matrix and a vector
  coAv = matrix(NA, nrow = nrow(expression), ncol = length(go))
  rownames(coAv) = rownames(expression)
  colnames(coAv) = names(go)

  auc = c()

  for (i in 1:steps) {

    exprStep = nEx[(i + size * (i - 1) ):(min(i + size * i, total)),]

    # making the correlation between all genes and genes in the text file
    stepCorrelation = cor(x = t(exprStep), y = t(libEx))
    rownames(stepCorrelation) = rownames(exprStep)
    colnames(stepCorrelation) = rownames(libEx)

    print(i)

    stepCorrelation[which(stepCorrelation == 1)] = NaN
    # removes correlation between protein/gene and itself
    # for (m in 1:length(co)) {

      # v1 = as.vector(arrayInd(m, dim(co)))

      # if (rownames(co)[v1[1]] == colnames(co)[v1[2]]) {
        # co[m] = NA
      # }

    # }

    # now calculate the means of each gene's correlation to each of the gene sets
    co_means = matrix(NA, nrow = nrow(stepCorrelation), ncol = length(go))
    rownames(co_means) = rownames(stepCorrelation)
    colnames(co_means) = names(go)

    # more efficient matrix making - taking a subset of genes based on gene set and 
    # averaging the correlations all at once
    for (j in 1:length(go)) {
      gene_set = go[[j]]
      
      if (length(gene_set) == 1) {
        co_means[, j] = stepCorrelation[, gene_set]
      } else {
        # now we subset and take the row means
        co_means[, j] = rowMeans(stepCorrelation[, gene_set], na.rm = TRUE)
      }
    }

    coAv[ (i+size * (i-1)) : (min(i+size * i, total)) , ] = co_means

    # so we construct the auc vector here - remember, the auc vector must be initialized outside of 
    # the giant step loop
    # now apply to all rows/genes
    for (current_gene in rownames(co_means)) {

      # this will eventually identify which sets have the gene and which sets do not
      containing_gene = co_means[current_gene,]
      
      # sorts the correlation averages
      containing_gene = sort(containing_gene, decreasing = TRUE)
      
      # for each gene_set in which the gene is a member, we want a 1 to be associated; 
      # for the rest, a zero.
      for (k in 1:length(go)) {
        if (current_gene %in% go[[k]]) {
          containing_gene[names(go)[k]] = 1
        } else {
          containing_gene[names(go)[k]] = 0
        }
      }
        
      cumulative = cumsum(containing_gene)
        
      # Part II: Now going to make a vector of AUC for all genes. 
      
      scaled_y = scale_vector(cumulative, 0, 1)
      scaled_x = scale_vector(1:length(cumulative), 0, 1)
      
      auc[current_gene] = trapz(scaled_x, scaled_y)
    }

  }

  # old go is necessary to revert to two-col form for Venn diagram creation
  ppiOld[[length(ppiOld) + 1]] = go
  names(ppiOld)[length(ppiOld)] = gmtNames[a]

  # make new GMT format out of correlation averages
  # will be put into two-col form at later point in another script
  toGmt_size(coAv, paste0(fileNames[a], "newGmt.txt"), go)

  # save the area under curve
  write.csv(auc, paste0("tempPpiAUC", fileNames[a], ".csv"))

  print(mean(auc, na.rm = T))

}
  