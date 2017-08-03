# one way to identify which libraries have these genes is to unlist the entire library and use which()
# this does not tell us which set the genes belong to. 
# indices = which(unlist(newLib) %in% genesOfInterest)

enrichrLibraries = c( # "Genes_Associated_with_NIH_Grants.txt", "Cancer_Cell_Line_Encyclopedia.txt", 
            "Achilles_fitness_decrease.txt", "Achilles_fitness_increase.txt", "Aging_Perturbations_from_GEO_down.txt", "Aging_Perturbations_from_GEO_up.txt", 
            # "Allen_Brain_Atlas_down.txt", "Allen_Brain_Atlas_up.txt", 
            "BioCarta_2013.txt", "BioCarta_2015.txt", "BioCarta_2016.txt", 
            # "BioPlex_2017.txt", 
            "ChEA_2013.txt", "ChEA_2015.txt", "ChEA_2016.txt", "Chromosome_Location.txt", 
            # "CORUM.txt", 
            "dbGaP.txt", "Disease_Perturbations_from_GEO_down.txt", "Disease_Perturbations_from_GEO_up.txt", 
            "Disease_Signatures_from_GEO_down_2014.txt", "Disease_Signatures_from_GEO_up_2014.txt", "Drug_Perturbations_from_GEO_2014.txt", 
            "Drug_Perturbations_from_GEO_down.txt", "Drug_Perturbations_from_GEO_up.txt", 
            # "DrugMatrix.txt", 
            "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt", "ENCODE_Histone_Modifications_2013.txt", "ENCODE_Histone_Modifications_2015.txt", 
            "ENCODE_TF_ChIP-seq_2014.txt", 
            "ENCODE_TF_2015.txt", 
            "Epigenomics_Roadmap_HM_ChIP-seq.txt", "ESCAPE.txt", 
            # "GeneSigDB.txt", 
            "Genome_Browser_PWMs.txt", "GO_Biological_Process_2013.txt", "GO_Biological_Process_2015.txt", 
            "GO_Biological_Process_2017.txt", "GO_Cellular_Component_2013.txt", "GO_Cellular_Component_2015.txt", 
            "GO_Cellular_Component_2017.txt", "GO_Molecular_Function_2013.txt", "GO_Molecular_Function_2015.txt", 
            "GO_Molecular_Function_2017.txt", 
            # "GTEx_Tissue_Sample_Gene_Expression_Profiles_down.txt", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up.txt", "HMDB_Metabolites.txt", 
            "HomoloGene.txt", "Human_Gene_Atlas.txt", 
            "Human_Phenotype_Ontology.txt", "HumanCyc_2015.txt", "Humancyc_2016.txt", "huMAP.txt", 
            # "Jensen_COMPARTMENTS.txt", "Jensen_DISEASES.txt", "Jensen_TISSUES.txt", 
            "KEA_2013.txt", "KEA_2015.txt", "KEGG_2013.txt", "KEGG_2015.txt", 
            "KEGG_2016.txt", "Kinase_Perturbations_from_GEO_down.txt", "Kinase_Perturbations_from_GEO_up.txt", 
            "Ligand_Perturbations_from_GEO_down.txt", "Ligand_Perturbations_from_GEO_up.txt", 
            # "LINCS_L1000_Chem_Pert_down.txt", "LINCS_L1000_Chem_Pert_up.txt", "LINCS_L1000_Kinase_Perturbations_down.txt", "LINCS_L1000_Kinase_Perturbations_up.txt", 
            "LINCS_L1000_Ligand_Perturbations_down.txt", "LINCS_L1000_Ligand_Perturbations_up.txt", 
            "MCF7_Perturbations_from_GEO_down.txt", "MCF7_Perturbations_from_GEO_up.txt", "MGI_Mammalian_Phenotype_2013.txt", 
            "MGI_Mammalian_Phenotype_2017.txt", "MGI_Mammalian_Phenotype_Level_3.txt", "MGI_Mammalian_Phenotype_Level_4.txt", 
            "Microbe_Perturbations_from_GEO_down.txt", "Microbe_Perturbations_from_GEO_up.txt", "Mouse_Gene_Atlas.txt", 
            "MSigDB_Computational.txt", "MSigDB_Oncogenic_Signatures.txt", "NCI-60_Cancer_Cell_Lines.txt", "NCI-Nature_2015.txt", 
            "NCI-Nature_2016.txt", "NURSA_Human_Endogenous_Complexome.txt", 
            # "Old_CMAP_down.txt", "Old_CMAP_up.txt", 
            "OMIM_Disease.txt", "OMIM_Expanded.txt", "Panther_2015.txt", "Panther_2016.txt", "Pfam_InterPro_Domains.txt", 
            "Phosphatase_Substrates_from_DEPOD.txt", "PPI_Hub_Proteins.txt", 
            "Reactome_2013.txt", "Reactome_2015.txt", "Reactome_2016.txt", 
            # "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO.txt", "SILAC_Phosphoproteomics.txt", 
            # "Single_Gene_Perturbations_from_GEO_down.txt", "Single_Gene_Perturbations_from_GEO_up.txt", 
            "TargetScan_microRNA.txt", "TF-LOF_Expression_from_GEO.txt", "Tissue_Protein_Expression_from_Human_Proteome_Map.txt", 
            "Tissue_Protein_Expression_from_ProteomicsDB.txt", "Transcription_Factor_PPIs.txt", "TRANSFAC_and_JASPAR_PWMs.txt", 
            "Virus_Perturbations_from_GEO_down.txt", "Virus_Perturbations_from_GEO_up.txt", "VirusMINT.txt", 
            "WikiPathways_2013.txt", 
            "WikiPathways_2015.txt", "WikiPathways_2016.txt")

# this is currently for gpcrs only
for (library in enrichrLibraries) {

	newLibrary = paste0("newGmt(3)_", library)
	oldLibrary = library 

	# analysis of new library
	# function makes text file GMT into list in R
	l.new = txtToList(newLibrary)
	# function gives us the genes from the genesOfInterest that are specifically in that library
	g.new = specify_genes(gpcr, l.new)
	# function gives us the sets that genes are in regarding a specific group of genes and a single library
	m.new = mark_sets(g.new, l.new)

	# analysis of old library
	# function makes text file GMT into list in R
	l.old = txtToList(oldLibrary)
	# function gives us the genes from the genesOfInterest that are specifically in that library
	g.old = specify_genes(gpcr, l.old)
	# function gives us the sets that genes are in regarding a specific group of genes and a single library
	m.old = mark_sets(g.old, l.old)

	# this function: it will tell us which sets contain the genes of interest in the new
	# BUT NOT in the old - this is for a specific library and a specific group of genes 
	# The current plan is that it returns a list containing all the genes in both the genes of interest and in the 
	# library, with each gene naming a vector of the sets that it is in.
	prediction1 = get_newOnly(gpcr, g.old, g.new, m.old, m.new)

	# check before using
	save(prediction1, file = paste0(gsub('.{4}$', '', library), "_predicted_gpcr_z.rda"))

}

for (gene in gpcr) {

      geneList = vector("list", length(enrichrLibraries))
      names(geneList) =  gsub('.{4}$', '', enrichrLibraries)

      for (library in enrichrLibraries) {

      	load(paste0(gsub('.{4}$', '', library), "_predicted_gpcr_z.rda"))

      	libName = gsub('.{4}$', '', library)

      	geneList[[libName]] = prediction1[[gene]]

      }

      # save(geneList, file = paste0(gene, "_predicted_z.rda"))
      sink(paste0("predicted_z_", gene, ".txt"))
      print(geneList)
      sink()
}