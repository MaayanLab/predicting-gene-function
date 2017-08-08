for (gene in ionChannel) {

      geneList = vector("list", length(enrichrLibraries))
      names(geneList) =  gsub('.{4}$', '', enrichrLibraries)

      for (library in enrichrLibraries) {

      	libName = gsub('.{4}$', '', library)

		load(paste0("predicted_z_ionChannel_rda/", libName, "_predicted_ionChannel_z.rda"))

      	geneList[[libName]] = prediction1[[gene]]

      }

      sink(paste0("gmtSite/predicted_z_", gene, ".txt"))

	for (a in 1:length(geneList)) {

		cat(names(geneList)[a], paste0(paste(geneList[[a]], collapse = "\t")), sep = "\t\t")
		cat("\n")
  
	}

	sink()
}