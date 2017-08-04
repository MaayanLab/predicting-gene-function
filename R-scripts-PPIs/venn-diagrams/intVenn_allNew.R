# # say we have goHu, goPlex, newGmtGrid - the old GMTs
# # newGmt signifies the new GMTs

# # Part II: trying to find common interactions between new GRID, new Plex, and new huMAP

# # Part I: Do common interactions between all 3.

# # gives us all the common protein sets between all three (old, new, new)
# toKeep = intersect(intersect(names(newGmtGrid), names(newGmtHu)), names(newGmtPlex))

# # gives us three GMTs in list format with the same protein set names in the same order
# newGmtGrid.abc = newGmtGrid[toKeep]
# newGmtHu.abc = newGmtHu[toKeep]
# newGmtPlex.abc = newGmtPlex[toKeep]

# intersection = c() # intersection of all three

# # onlyGrid = c() # only in the old BioGRID
# # onlyHu = c() # only in the new huMAP
# # onlyPlex = c() # only in the new BioPlex

# for (u in 1:length(toKeep)) {

# 	intersection[u] = length(intersect(intersect(newGmtGrid.abc[[u]], newGmtHu.abc[[u]]), newGmtPlex.abc[[u]]))

# 	# onlyGrid[u] = length(setdiff(setdiff(newGmtGrid.s[[u]], newGmtHu.s[[u]]), newGmtPlex.s[[u]]))
# 	# onlyHu[u] = length(setdiff(setdiff(newGmtHu.s[[u]], newGmtGrid.s[[u]]), newGmtPlex.s[[u]]))
# 	# onlyPlex[u] = length(setdiff(setdiff(newGmtPlex.s[[u]], newGmtGrid.s[[u]]), newGmtHu.s[[u]]))

# 	# GridHu[u] = length(setdiff(intersect(newGmtGrid.s[[u]], newGmtHu.s[[u]]), newGmtPlex.s[[u]]))
# 	# HuPlex[u] = length(setdiff(intersect(newGmtPlex.s[[u]], newGmtHu.s[[u]]), newGmtGrid.s[[u]]))
# 	# GridPlex[u] = length(setdiff(intersect(newGmtGrid.s[[u]], newGmtPlex.s[[u]]), newGmtHu.s[[u]]))

# }

# # Part II: Do common interactions between pairs of 2, including the interactions between all 3.

# GridHu = c() # both new BioGRID and new huMAP (including all 3)
# HuPlex = c() # both new BioPlex and new huMAP
# GridPlex = c() # both new BioGRID and new BioPlex

# # between A and B
# toKeep = intersect(names(newGmtGrid), names(newGmtHu))

# newGmtGrid.ab = newGmtGrid[toKeep]
# newGmtHu.ab = newGmtHu[toKeep]

# for (u in 1:length(toKeep)) {

# 	GridHu[u] = length(intersect(newGmtGrid.ab[[u]], newGmtHu.ab[[u]]))
# }

# # between B and C
# toKeep = intersect(names(newGmtHu), names(newGmtPlex))

# newGmtHu.bc = newGmtHu[toKeep]
# newGmtPlex.bc = newGmtPlex[toKeep]

# for (u in 1:length(toKeep)) {

# 	HuPlex[u] = length(intersect(newGmtHu.bc[[u]], newGmtPlex.bc[[u]]))
# }

# # between A and C
# toKeep = intersect(names(newGmtGrid), names(newGmtPlex))

# newGmtGrid.ac = newGmtGrid[toKeep]
# newGmtPlex.ac = newGmtPlex[toKeep]

# for(u in 1:length(toKeep)) {

# 	GridPlex[u] = length(intersect(newGmtPlex.ac[[u]], newGmtGrid.ac[[u]]))
# }


# sumIntersect = sum(intersection)

# # sumOnlyGrid = sum(onlyGrid)
# # sumOnlyHu = sum(onlyHu)
# # sumOnlyPlex = sum(onlyPlex)

# sumGridHu = sum(GridHu)
# sumHuPlex = sum(HuPlex)
# sumGridPlex = sum(GridPlex)

# draw.triple.venn(sum(lengths(newGmtGrid)), 
# 	             sum(lengths(newGmtHu)),
# 	             sum(lengths(newGmtPlex)),
# 	             sumGridHu,
# 	             sumHuPlex,
# 	             sumGridPlex,
# 	             sumIntersect,
# 	             category = c("New GMT, BioGRID", "New GMT, huMAP", "New GMT, BioPlex"),
# 			 scaled = FALSE,
# 			 fill = c("goldenrod3", "dodgerblue", "chartreuse3"),
# 			 lty = "blank",
# 			 cex = rep(3, 7),
# 			 cat.pos = c(-15, 15, 180),
# 			 cat.dist = c(0.04, 0.04, 0.03),
# 			 cat.cex = c(2, 2, 2)
# 			 )

require(data.table)

nG = setkey(data.table(newGrid))
nH = setkey(data.table(newHu))
nG.nH = length(na.omit(nH[nG, which=TRUE]))

nH = setkey(data.table(newHu))
nP = setkey(data.table(newPlex))
nH.nP = length(na.omit(nP[nH, which=TRUE]))

nG = setkey(data.table(newGrid))
nP = setkey(data.table(newPlex))
nG.nP = length(na.omit(nP[nG, which=TRUE]))

int1 = na.omit(nP[nH, which=TRUE])

int2 = length(na.omit(nG[nP[int1, ], which = TRUE]))

draw.triple.venn(nrow(newGrid), 
	nrow(newHu),
	nrow(newPlex),
	nG.nH,
	nH.nP,
	nG.nP,
	int2,
	category = c("New PPI, BioGRID", "New PPI, huMAP", "New PPI, BioPlex"),
	scaled = FALSE,
	fill = c("darkgoldenrod1", "dodgerblue", "chartreuse3"),
	lty = "blank",
	cex = rep(3, 7),
	cat.pos = c(-15, 15, 180),
	cat.dist = c(0.04, 0.04, 0.03),
	cat.cex = c(2, 2, 2)
	)

