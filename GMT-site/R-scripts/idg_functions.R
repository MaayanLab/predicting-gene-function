kinase = unlist(strsplit(readLines("idg_u24_kinases.txt"), ", "))
ionChannel = unlist(strsplit(readLines("idg_u24_ionChan.txt"), ", "))
gpcr = unlist(strsplit(readLines("idg_u24_gpcr.txt"), ", "))

# one way to identify which libraries have these genes is to unlist the entire library and use which()
# this does not tell us which set the genes belong to. 
# indices = which(unlist(newLib) %in% genesOfInterest)

# function makes text file GMT into list in R
txtToList <- function(txtFile) {
  gmt = readLines(txtFile)
  
  go = list()
  for(line in gmt){
    # splits line into the tabbed words
    sp = unlist(strsplit(line, "\t"))
    
    # term is the name of the current set
    term = sp[1]
    
    # t ends up being all the current set's genes, but it's currently empty
    t = c()
    
    # starts at 3 because that's where the genes begin (the name of the set, an empty space, then first gene)
    for(m in 3:length(sp)){
      sp1 = unlist(strsplit( sp[m], ","))
      t = c(t, sp1[1])
    }
    # holds all the phenotypes and their genes in the list format, where the name of the set is the category name
    # Each set is a vector of its genes
    # taking the intersection of the names in expression and the names in this particular list element/set
    go[[length(go)+1]] = t 
    names(go)[length(go)] = term
  }
  
  return(go)
}

# function gives us the genes from the genesOfInterest that are specifically in that library
specify_genes <- function(genesOfInterest, lib) {
	g1 = genesOfInterest[which(genesOfInterest %in% unlist(lib))]
	# print(paste0("Library contains ", length(g1), " genes out of the genes of interest."))
	return(g1)
}

# function gives us the sets that genes are in regarding a specific group of genes and a single library
mark_sets <- function(g1, lib) {

	m1 = matrix(NA, nrow = length(g1), ncol = length(lib), dimnames = list(g1, names(lib)))

	# for each set in the new library, find which genes from the pertinent genes of interest (g1) are in each set
	# take the genes that are in the set a and mark them as "TRUE"
	for (a in 1:length(lib)) {
		v1 = which(g1 %in% lib[[a]])
		m1[v1, a] = T
	}

	return(m1)
}

# we'll have newBin as the binary matrix of the new library and oldBin as the binary matrix of the old library


# The eventual hope for this function is that it will tell us which sets contain the genes of interest in the new
# BUT NOT in the old - this is for a specific library and a specific group of genes 
# The current plan is that it returns a list containing all the genes in both the genes of interest and in the 
# library, with each gene naming a vector of the sets that it is in.
get_newOnly <- function(focus, g1.old, g1.new, oldBin, newBin) {

	g2 = c(intersect(g1.old, g1.new), g1.new)
	g2 = unique(g2)
	
	predictions = vector("list", length(focus))
	names(predictions) = focus

	for (current_gene in focus) {

		if (current_gene %in% g2) {
			# newSets - for the current_gene, these are the sets that the gene is in regarding the new library
			newSets = names(which(newBin[current_gene, ]))

			if (current_gene %in% g1.old) {
				# oldSets - for the current_gene, these are the sets that the gene is in regarding the old library 
				oldSets = names(which(oldBin[current_gene, ]))

				if (length(setdiff(newSets, oldSets)) != 0) {
					# get the sets that are in new but not old
					predictions[[current_gene]] = c(predictions[[current_gene]], setdiff(newSets, oldSets))
				} else {
					predictions[[current_gene]] = "No new information predicted"
				}
			} else {
				predictions[[current_gene]] = newSets
			}
		} else {
			predictions[[current_gene]] = "Gene is not in new version of this library"
		}
	}

	return(predictions)
}
