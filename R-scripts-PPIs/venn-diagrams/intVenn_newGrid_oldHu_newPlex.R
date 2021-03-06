# # say we have goHu, goPlex, newGmtGrid - the old GMTs
# # newGmt signifies the new GMTs

# # Part II: trying to find common interactions between new GRID, old HuMAP, and new Plex

# # Part I: Do common interactions between all 3.

# # gives us all the common protein sets between all three (new, old, new)
# toKeep = intersect(intersect(names(newGmtGrid), names(goHu)), names(newGmtPlex))

# # gives us three GMTs in list format with the same protein set names in the same order
# newGmtGrid.abc = newGmtGrid[toKeep]
# goHu.abc = goHu[toKeep]
# newGmtPlex.abc = newGmtPlex[toKeep]

# intersection = c() # intersection of all three

# for (u in 1:length(toKeep)) {

# 	intersection[u] = length(intersect(intersect(newGmtGrid.abc[[u]], goHu.abc[[u]]), newGmtPlex.abc[[u]]))

# }

# # Part II: Do common interactions between pairs of 2, including the interactions between all 3.

# GridHu = c() # both new BioGRID and old huMAP (including all 3)
# HuPlex = c() # both new BioPlex and old huMAP
# GridPlex = c() # both new BioGRID and new BioPlex

# # between A and B
# toKeep = intersect(names(newGmtGrid), names(goHu))

# newGmtGrid.ab = newGmtGrid[toKeep]
# goHu.ab = goHu[toKeep]

# for (u in 1:length(toKeep)) {

# 	GridHu[u] = length(intersect(newGmtGrid.ab[[u]], goHu.ab[[u]]))
# }

# # between B and C
# toKeep = intersect(names(goHu), names(newGmtPlex))

# goHu.bc = goHu[toKeep]
# newGmtPlex.bc = newGmtPlex[toKeep]

# for (u in 1:length(toKeep)) {

# 	HuPlex[u] = length(intersect(goHu.bc[[u]], newGmtPlex.bc[[u]]))
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
# 	             sum(lengths(goHu)),
# 	             sum(lengths(newGmtPlex)),
# 	             sumGridHu,
# 	             sumHuPlex,
# 	             sumGridPlex,
# 	             sumIntersect,
# 	             category = c("New GMT, BioGRID", "Original GMT, huMAP", "New GMT, BioPlex"),
# 			 scaled = FALSE,
# 			 fill = c("goldenrod3", "lightskyblue", "chartreuse3"),
# 			 lty = "blank",
# 			 cex = rep(3, 7),
# 			 cat.pos = c(-15, 15, 180),
# 			 cat.dist = c(0.04, 0.04, 0.03),
# 			 cat.cex = c(2, 2, 2)
# 			 )

require(data.table)

nG = setkey(data.table(newGrid))
oH = setkey(data.table(oldHu))
nG.oH = length(na.omit(oH[nG, which=TRUE]))

oH = setkey(data.table(oldHu))
nP = setkey(data.table(newPlex))
oH.nP = length(na.omit(nP[oH, which=TRUE]))

nG = setkey(data.table(newGrid))
nP = setkey(data.table(newPlex))
nG.nP = length(na.omit(nP[nG, which=TRUE]))

inter = na.omit(nP[oH, which=TRUE])

intersection = length(na.omit(nG[nP[inter, ], which = TRUE]))

draw.triple.venn(nrow(newGrid), 
	nrow(oldHu),
	nrow(newPlex),
	nG.oH,
	oH.nP,
	nG.nP,
	intersection,
	category = c("New PPI, BioGRID", "Original PPI, huMAP", "New PPI, BioPlex"),
	scaled = FALSE,
	fill = c("darkgoldenrod1", "lightskyblue", "chartreuse3"),
	lty = "blank",
	cex = rep(3, 7),
	cat.pos = c(-15, 15, 180),
	cat.dist = c(0.04, 0.04, 0.03),
	cat.cex = c(2, 2, 2)
	)
