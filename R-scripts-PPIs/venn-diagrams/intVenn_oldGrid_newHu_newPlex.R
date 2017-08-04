# # say we have goHu, goPlex, goGrid - the old GMTs
# # newGmt signifies the new GMTs

# # Part II: trying to find common interactions between old GRID, new Plex, and new huMAP

# # toKeep3 = intersect(names(goHu), names(goPlex))
# # toKeep3 = intersect(toKeep3, names(goGrid))

# # goHu3 = goHu[toKeep3]
# # goPlex3 = goPlex[toKeep3]
# # goGrid3 = goGrid[toKeep3]

# # Part I: Do common interactions between all 3.

# # gives us all the common protein sets between all three (old, new, new)
# toKeep = intersect(intersect(names(goGrid), names(newGmtHu)), names(newGmtPlex))

# # gives us three GMTs in list format with the same protein set names in the same order
# goGrid.abc = goGrid[toKeep]
# newGmtHu.abc = newGmtHu[toKeep]
# newGmtPlex.abc = newGmtPlex[toKeep]

# intersection = c() # intersection of all three

# # onlyGrid = c() # only in the old BioGRID
# # onlyHu = c() # only in the new huMAP
# # onlyPlex = c() # only in the new BioPlex

# for (u in 1:length(toKeep)) {

# 	intersection[u] = length(intersect(intersect(goGrid.abc[[u]], newGmtHu.abc[[u]]), newGmtPlex.abc[[u]]))

# 	# onlyGrid[u] = length(setdiff(setdiff(goGrid.s[[u]], newGmtHu.s[[u]]), newGmtPlex.s[[u]]))
# 	# onlyHu[u] = length(setdiff(setdiff(newGmtHu.s[[u]], goGrid.s[[u]]), newGmtPlex.s[[u]]))
# 	# onlyPlex[u] = length(setdiff(setdiff(newGmtPlex.s[[u]], goGrid.s[[u]]), newGmtHu.s[[u]]))

# 	# GridHu[u] = length(setdiff(intersect(goGrid.s[[u]], newGmtHu.s[[u]]), newGmtPlex.s[[u]]))
# 	# HuPlex[u] = length(setdiff(intersect(newGmtPlex.s[[u]], newGmtHu.s[[u]]), goGrid.s[[u]]))
# 	# GridPlex[u] = length(setdiff(intersect(goGrid.s[[u]], newGmtPlex.s[[u]]), newGmtHu.s[[u]]))

# }

# # Part II: Do common interactions between pairs of 2, including the interactions between all 3.

# GridHu = c() # both old BioGRID and new huMAP (including all 3)
# HuPlex = c() # both new BioPlex and new huMAP
# GridPlex = c() # both old BioGRID and new BioPlex

# # between A and B
# toKeep = intersect(names(goGrid), names(newGmtHu))

# goGrid.ab = goGrid[toKeep]
# newGmtHu.ab = newGmtHu[toKeep]

# for (u in 1:length(toKeep)) {

# 	GridHu[u] = length(intersect(goGrid.ab[[u]], newGmtHu.ab[[u]]))
# }

# # between B and C
# toKeep = intersect(names(newGmtHu), names(newGmtPlex))

# newGmtHu.bc = newGmtHu[toKeep]
# newGmtPlex.bc = newGmtPlex[toKeep]

# for (u in 1:length(toKeep)) {

# 	HuPlex[u] = length(intersect(newGmtHu.bc[[u]], newGmtPlex.bc[[u]]))
# }

# # between A and C
# toKeep = intersect(names(goGrid), names(newGmtPlex))

# goGrid.ac = goGrid[toKeep]
# newGmtPlex.ac = newGmtPlex[toKeep]

# for(u in 1:length(toKeep)) {

# 	GridPlex[u] = length(intersect(newGmtPlex.ac[[u]], goGrid.ac[[u]]))
# }


# sumIntersect = sum(intersection)

# # sumOnlyGrid = sum(onlyGrid)
# # sumOnlyHu = sum(onlyHu)
# # sumOnlyPlex = sum(onlyPlex)

# sumGridHu = sum(GridHu)
# sumHuPlex = sum(HuPlex)
# sumGridPlex = sum(GridPlex)

require(data.table)

oG = setkey(data.table(oldGrid))
nH = setkey(data.table(newHu))
oG.nH = length(na.omit(nH[oG, which=TRUE]))

nH = setkey(data.table(newHu))
nP = setkey(data.table(newPlex))
nH.nP = length(na.omit(nP[nH, which=TRUE]))

oG = setkey(data.table(oldGrid))
nP = setkey(data.table(newPlex))
oG.nP = length(na.omit(nP[oG, which=TRUE]))

inter = na.omit(nP[nH, which=TRUE])

intersection = length(na.omit(oG[nP[inter, ], which = TRUE]))

draw.triple.venn(nrow(oldGrid), 
	nrow(newHu),
	nrow(newPlex),
	oG.nH,
	nH.nP,
	oG.nP,
	intersection,
	category = c("Original PPI, BioGRID", "New PPI, huMAP", "New PPI, BioPlex"),
	scaled = FALSE,
	fill = c("gold1", "dodgerblue", "chartreuse3"),
	lty = "blank",
	cex = rep(3, 7),
	cat.pos = c(-15, 15, 180),
	cat.dist = c(0.04, 0.04, 0.03),
	cat.cex = c(2, 2, 2)
	)

