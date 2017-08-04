# # say we have newGmtHu, goPlex, newGmtGrid - the old GMTs
# # newGmt signifies the new GMTs

# # Part II: trying to find common interactions between new GRID, new HuMAP, and old Plex

# # Part I: Do common interactions between all 3.

# # gives us all the common protein sets between all three (new, new, old)
# toKeep = intersect(intersect(names(newGmtGrid), names(newGmtHu)), names(goPlex))

# # gives us three GMTs in list format with the same protein set names in the same order
# newGmtGrid.abc = newGmtGrid[toKeep]
# newGmtHu.abc = newGmtHu[toKeep]
# goPlex.abc = goPlex[toKeep]

# intersection = c() # intersection of all three

# for (u in 1:length(toKeep)) {

# 	intersection[u] = length(intersect(intersect(newGmtGrid.abc[[u]], newGmtHu.abc[[u]]), goPlex.abc[[u]]))

# }

# # Part II: Do common interactions between pairs of 2, including the interactions between all 3.

# GridHu = c() # both new BioGRID and old huMAP (including all 3)
# HuPlex = c() # both new BioPlex and old huMAP
# GridPlex = c() # both new BioGRID and new BioPlex

# # between A and B
# toKeep = intersect(names(newGmtGrid), names(newGmtHu))

# newGmtGrid.ab = newGmtGrid[toKeep]
# newGmtHu.ab = newGmtHu[toKeep]

# for (u in 1:length(toKeep)) {

# 	GridHu[u] = length(intersect(newGmtGrid.ab[[u]], newGmtHu.ab[[u]]))
# }

# # between B and C
# toKeep = intersect(names(newGmtHu), names(goPlex))

# newGmtHu.bc = newGmtHu[toKeep]
# goPlex.bc = goPlex[toKeep]

# for (u in 1:length(toKeep)) {

# 	HuPlex[u] = length(intersect(newGmtHu.bc[[u]], goPlex.bc[[u]]))
# }

# # between A and C
# toKeep = intersect(names(newGmtGrid), names(goPlex))

# newGmtGrid.ac = newGmtGrid[toKeep]
# goPlex.ac = goPlex[toKeep]

# for(u in 1:length(toKeep)) {

# 	GridPlex[u] = length(intersect(goPlex.ac[[u]], newGmtGrid.ac[[u]]))
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
# 	             sum(lengths(goPlex)),
# 	             sumGridHu,
# 	             sumHuPlex,
# 	             sumGridPlex,
# 	             sumIntersect,
# 	             category = c("New GMT, BioGRID", "New GMT, huMAP", "Original GMT, BioPlex"),
# 			 scaled = FALSE,
# 			 fill = c("darkgoldenrod1", "dodgerblue", "lawngreen"),
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
oP = setkey(data.table(oldPlex))
nH.oP = length(na.omit(oP[nH, which=TRUE]))

nG = setkey(data.table(newGrid))
oP = setkey(data.table(oldPlex))
nG.oP = length(na.omit(oP[nG, which=TRUE]))

inter = na.omit(oP[nH, which=TRUE])

intersection = length(na.omit(nG[oP[inter, ], which = TRUE]))

draw.triple.venn(nrow(newGrid), 
	nrow(newHu),
	nrow(oldPlex),
	nG.nH,
	nH.oP,
	nG.oP,
	intersection,
	category = c("New PPI, BioGRID", "New PPI, huMAP", "Original PPI, BioPlex"),
	scaled = FALSE,
	fill = c("darkgoldenrod1", "dodgerblue", "lawngreen"),
	lty = "blank",
	cex = rep(3, 7),
	cat.pos = c(-15, 15, 180),
	cat.dist = c(0.04, 0.04, 0.03),
	cat.cex = c(2, 2, 2)
	)

