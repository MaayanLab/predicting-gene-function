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

destination_file = "CCLE_Expression_Entrez_2012-09-29.gct"
extracted_expression_file = "sample_expression_matrix.tsv"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
  print("Downloading compressed gene expression matrix.")
  url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
  download(url, destination_file, mode="wb")
} else{
  print("Local file already exists.")
}

gct = readLines("CCLE_Expression_Entrez_2012-09-29.gct")

sp = unlist(strsplit(gct[2], "\t"))
ccle = matrix(NA, as.numeric(sp[1]), as.numeric(sp[2]))

sp = unlist(strsplit(gct[3], "\t"))
colnames(ccle) = sp[3:length(sp)]

genes = c()

for(n in 4:length(gct)) {

  # splits line into the tabbed words
  sp = unlist(strsplit(gct[n], "\t"))
  genes[n - 3] = sp[2] 

  ccle[n - 3, ] = as.numeric(sp[3:length(sp)]) 

}

rownames(ccle) = genes

toRm = which(rownames(ccle) %in% "")
ccle = ccle[-toRm, ]

# the matrix version of expression creation
expression = ccle
# genes
rownames(expression) = rownames(ccle)
# samples
colnames(expression) = colnames(ccle)


# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))

cat("Dimensions of Expression: ", dim(expression), "\n")

# we want to do quantile normalization on the entire expression matrix
nEx = normalize.quantiles(expression)
rownames(nEx) = rownames(ccle)
colnames(nEx) = colnames(ccle)

# Part II: Make the expression matrix only the genes also in the gmt file. To do so, read the gmt file.

gmtNames = c("ChEA_2016.txt", "ENCODE_TF_2015.txt", "GO_Biological_Process_2017.txt", "GO_Cellular_Component_2017.txt", "GO_Molecular_Function_2017.txt", "KEA_2015.txt", "KEGG_2016.txt", "MGI_Mammalian_Phenotype_Level_4.txt")

# gmtNames = c("Genes_Associated_with_NIH_Grants.txt", "Cancer_Cell_Line_Encyclopedia.txt", "Achilles_fitness_decrease.txt", "Achilles_fitness_increase.txt", "Aging_Perturbations_from_GEO_down.txt", "Aging_Perturbations_from_GEO_up.txt",
#   "Allen_Brain_Atlas_down.txt", "Allen_Brain_Atlas_up.txt", "BioCarta_2013.txt", "BioCarta_2015.txt", "BioCarta_2016.txt", "BioPlex_2017.txt", "ChEA_2013.txt", "ChEA_2015.txt", "ChEA_2016.txt", "Chromosome_Location.txt", "CORUM.txt", 
#   "dbGaP.txt", "Disease_Perturbations_from_GEO_down.txt", "Disease_Perturbations_from_GEO_up.txt", "Disease_Signatures_from_GEO_down_2014.txt", "Disease_Signatures_from_GEO_up_2014.txt", "Drug_Perturbations_from_GEO_2014.txt", 
#   "Drug_Perturbations_from_GEO_down.txt", "Drug_Perturbations_from_GEO_up.txt", "DrugMatrix.txt", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt", "ENCODE_Histone_Modifications_2013.txt", "ENCODE_Histone_Modifications_2015.txt", 
#   "ENCODE_TF_ChIP-seq_2014.txt", "ENCODE_TF_ChIP-seq_2015.txt", "Epigenomics_Roadmap_HM_ChIP-seq.txt", "ESCAPE.txt", "GeneSigDB.txt", "Genome_Browser_PWMs.txt", "GO_Biological_Process_2013.txt", "GO_Biological_Process_2015.txt", 
#   "GO_Biological_Process_2017.txt", "GO_Cellular_Component_2013.txt", "GO_Cellular_Component_2015.txt", "GO_Cellular_Component_2017.txt", "GO_Molecular_Function_2013.txt", "GO_Molecular_Function_2015.txt", 
#   "GO_Molecular_Function_2017.txt", "GTEx_Tissue_Sample_Gene_Expression_Profiles_down.txt", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up.txt", "HMDB_Metabolites.txt", "HomoloGene.txt", "Human_Gene_Atlas.txt", 
#   "Human_Phenotype_Ontology.txt", "HumanCyc_2015.txt", "Humancyc_2016.txt", "huMAP.txt", "Jensen_COMPARTMENTS.txt", "Jensen_DISEASES.txt", "Jensen_TISSUES.txt", "KEA_2013.txt", "KEA_2015.txt", "KEGG_2013.txt", "KEGG_2015.txt", 
#   "KEGG_2016.txt", "Kinase_Perturbations_from_GEO_down.txt", "Kinase_Perturbations_from_GEO_up.txt", "Ligand_Perturbations_from_GEO_down.txt", "Ligand_Perturbations_from_GEO_up.txt", "LINCS_L1000_Chem_Pert_down.txt", 
#   "LINCS_L1000_Chem_Pert_up.txt", "LINCS_L1000_Kinase_Perturbations_down.txt", "LINCS_L1000_Kinase_Perturbations_up.txt", "LINCS_L1000_Ligand_Perturbations_down.txt", "LINCS_L1000_Ligand_Perturbations_up.txt", 
#   "MCF7_Perturbations_from_GEO_down.txt", "MCF7_Perturbations_from_GEO_up.txt", "MGI_Mammalian_Phenotype_2013.txt", "MGI_Mammalian_Phenotype_2017.txt", "MGI_Mammalian_Phenotype_Level_3.txt", "MGI_Mammalian_Phenotype_Level_4.txt", 
#   "Microbe_Perturbations_from_GEO_down.txt", "Microbe_Perturbations_from_GEO_up.txt", "Mouse_Gene_Atlas.txt", "MSigDB_Computational.txt", "MSigDB_Oncogenic_Signatures.txt", "NCI-60_Cancer_Cell_Lines.txt", "NCI-Nature_2015.txt", 
#   "NCI-Nature_2016.txt", "NURSA_Human_Endogenous_Complexome.txt", "Old_CMAP_down.txt", "Old_CMAP_up.txt", "OMIM_Disease.txt", "OMIM_Expanded.txt", "Panther_2015.txt", "Panther_2016.txt", "Pfam_InterPro_Domains.txt", 
#   "Phosphatase_Substrates_from_DEPOD.txt", "PPI_Hub_Proteins.txt", "Reactome_2013.txt", "Reactome_2015.txt", "Reactome_2016.txt", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO.txt", "SILAC_Phosphoproteomics.txt", 
#   "Single_Gene_Perturbations_from_GEO_down.txt", "Single_Gene_Perturbations_from_GEO_up.txt", "TargetScan_microRNA.txt", "TF-LOF_Expression_from_GEO.txt", "Tissue_Protein_Expression_from_Human_Proteome_Map.txt", 
#   "Tissue_Protein_Expression_from_ProteomicsDB.txt", "Transcription_Factor_PPIs.txt", "TRANSFAC_and_JASPAR_PWMs.txt", "Virus_Perturbations_from_GEO_down.txt", "Virus_Perturbations_from_GEO_up.txt", "VirusMINT, WikiPathways_2013.txt", 
#   "WikiPathways_2015.txt", "WikiPathways_2016.txt")


# This script just extracts all the genes from the H_P_O file and removes duplicates. It puts them into a vector.
# it first reads the gmt file

ccleCoAvg = list()
ccleAUC = list()

# HERE BEGINS THE FOR LOOP THAT PRODUCES AUC AND average correlation matrices FOR EACH OF THE GMT FILES 
for (txtFile in gmtNames) {

  print(txtFile)

  gmt = readLines(txtFile)

  go = list()
  for(line in gmt) {
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

  # a matrix and a vector
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
    # for (m in 1:length(stepCorrelation)) {

    #   v1 = as.vector(arrayInd(m, dim(stepCorrelation)))

    #   if (rownames(stepCorrelation)[v1[1]] == colnames(stepCorrelation)[v1[2]]) {
    #     stepCorrelation[m] = NA
    #   }

    # }
    
    # now calculate the means of each gene's correlation to each of the gene sets
    co_means = matrix(NA, nrow = nrow(stepCorrelation), ncol = length(names(go)))
    rownames(co_means) = rownames(stepCorrelation)
    colnames(co_means) = names(go)
    
    # more efficient matrix making - taking a subset of genes based on gene set and averaging the correlations all at once
    for (j in 1:length(go)) {
      gene_set = go[[j]]
      
      if (length(gene_set) == 1) {
        co_means[, j] = stepCorrelation[, gene_set]
      } else {
        co_means[, j] = rowMeans(stepCorrelation[,gene_set], na.rm = TRUE)
      }
    }
    
    allCrrAvg[ (i+size * (i-1)) : (min(i+size * i, total)) , ] = co_means
    
    # so we construct the auc vector here - remember, the auc vector must be initialized outside of the giant step loop
    # now apply to all rows/genes
    for (current_gene in rownames(co_means)) {
      
      # this will eventually identify which sets have the gene and which sets do not
      containing_gene = co_means[current_gene,]
      
      # sorts the correlation averages
      containing_gene = sort(containing_gene, decreasing = TRUE)
      
      # for each gene_set in which the gene is a member, we want a 1 to be associated; for the rest, a zero.
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

  ccleCoAvg[[length(ccleCoAvg) + 1]] = allCrrAvg
  names(ccleCoAvg)[length(ccleCoAvg)] = txtFile
  
  ccleAUC[[length(ccleAUC) + 1]] = auc
  names(ccleAUC)[length(ccleAUC)] = txtFile

  print(mean(auc, na.rm = T))

}