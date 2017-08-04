toGmt_size <- function(corrAvg, fileName, oldGmt) {

  newGmt = list()

  for (b in 1:ncol(corrAvg)) {
    setOne = (corrAvg)[,b] 
    names(setOne) = rownames(corrAvg)

    setOne = sort(setOne, decreasing = TRUE)

    newSet = setOne[1:length(oldGmt[[b]])]

    newGmt[[length(newGmt) + 1]] = unique(names(newSet))
    names(newGmt)[length(newGmt)] = colnames(corrAvg)[b]

  }
  sink(fileName)
  for (i in 1:length(newGmt)) {
    cat(names(newGmt)[i], paste0(paste(newGmt[[i]], collapse = "\t")), sep = "\t\t")
    cat("\n")
  }
  sink()
}


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


cat("Dimensions of Expression: ", dim(expression), "\n")

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

# HERE BEGINS THE FOR LOOP THAT PRODUCES AUC AND average correlation matrices FOR EACH OF THE GMT FILES 
for (a in 1:length(gmtNames)) {

  print(gmtNames[a])

  gmt = readLines(gmtNames[a])

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

  toGmt_size(ccleCoAvg[[a]], paste0("newCcle", gmtNames[a]), go)

}
