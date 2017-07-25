# newGmt <- list()

# for (n in 1:ncol(listCoAvg[[4]])) {
# 	setOne <- (listCoAvg[[4]])[,n] 
# 	names(setOne) <- rownames(listCoAvg[[4]])

# 	setOne <- sort(setOne)

# 	pop_mean <- mean(setOne)
# 	pop_sd <- sd(setOne)

# 	pop_z <- (setOne - pop_mean) / pop_sd

# 	# cutting off at 3 standard deviations
# 	# cutOff <- sort(pop_z[pop_z >= 3], decreasing = TRUE)
# 	newSet <- names(pop_z[pop_z >= 4]) 

# 	newGmt[[length(newGmt) + 1]] <- newSet
# 	names(newGmt)[length(newGmt)] <- colnames(listCoAvg[[4]])[n]

# }

# by standard deviation cut-off

toGmt_z <- function(corrAvg, fileName, z) {

	newGmt = list()

	for (n in 1:ncol(corrAvg)) {
		setOne = (corrAvg)[,n] 
		names(setOne) = rownames(corrAvg)

		setOne = sort(setOne, decreasing = TRUE)

		pop_mean = mean(setOne)
		pop_sd = sd(setOne)

		pop_z = (setOne - pop_mean) / pop_sd

	# cutting off at 3 standard deviations
	newSet = names(pop_z[pop_z >= z]) 

	newGmt[[length(newGmt) + 1]] = unique(newSet)
	names(newGmt)[length(newGmt)] = colnames(corrAvg)[n]

	}

# toRm = c()

# for (i in 1:length(newGmt)) {

	# if (length(newGmt[[i]]) == 0) {
		# toRm = c(toRm, i)
	# }

# }

# newGmt = newGmt[-toRm]

	sink(fileName)

	for (i in 1:length(newGmt)) {

		cat(names(newGmt)[i], paste0(paste(newGmt[[i]], collapse = "\t")), sep = "\t\t")
		cat("\n")

	}

	sink()

}
