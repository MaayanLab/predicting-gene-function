# # say we have goHu, goPlex, goGrid - the old GMTs
# # newGmt signifies the new GMTs

# # Part II: trying to find common interactions between old GRID, old Plex, and old huMAP

# # toKeep3 = intersect(names(goHu), names(goPlex))
# # toKeep3 = intersect(toKeep3, names(goGrid))

# # goHu3 = goHu[toKeep3]
# # goPlex3 = goPlex[toKeep3]
# # goGrid3 = goGrid[toKeep3]

# # Part I: Do common interactions between all 3.

# # gives us all the common protein sets between all three (old, new, new)
# toKeep = intersect(intersect(names(goGrid), names(goHu)), names(goPlex))

# # gives us three GMTs in list format with the same protein set names in the same order
# goGrid.abc = goGrid[toKeep]
# goHu.abc = goHu[toKeep]
# goPlex.abc = goPlex[toKeep]

# intersection = c() # intersection of all three

# for (u in 1:length(toKeep)) {

# 	intersection[u] = length(intersect(intersect(goGrid.abc[[u]], goHu.abc[[u]]), goPlex.abc[[u]]))

# 	# onlyGrid[u] = length(setdiff(setdiff(goGrid.s[[u]], goHu.s[[u]]), goPlex.s[[u]]))
# 	# onlyHu[u] = length(setdiff(setdiff(goHu.s[[u]], goGrid.s[[u]]), goPlex.s[[u]]))
# 	# onlyPlex[u] = length(setdiff(setdiff(goPlex.s[[u]], goGrid.s[[u]]), goHu.s[[u]]))

# 	# GridHu[u] = length(setdiff(intersect(goGrid.s[[u]], goHu.s[[u]]), goPlex.s[[u]]))
# 	# HuPlex[u] = length(setdiff(intersect(goPlex.s[[u]], goHu.s[[u]]), goGrid.s[[u]]))
# 	# GridPlex[u] = length(setdiff(intersect(goGrid.s[[u]], goPlex.s[[u]]), goHu.s[[u]]))

# }

# # Part II: Do common interactions between pairs of 2, including the interactions between all 3.

# GridHu = c() # both old BioGRID and old huMAP (including all 3)
# HuPlex = c() # both old BioPlex and old huMAP
# GridPlex = c() # both old BioGRID and old BioPlex

# # between A and B
# toKeep = intersect(names(goGrid), names(goHu))

# goGrid.ab = goGrid[toKeep]
# goHu.ab = goHu[toKeep]

# for (u in 1:length(toKeep)) {

# 	GridHu[u] = length(intersect(goGrid.ab[[u]], goHu.ab[[u]]))
# }

# # between B and C
# toKeep = intersect(names(goHu), names(goPlex))

# goHu.bc = goHu[toKeep]
# goPlex.bc = goPlex[toKeep]

# for (u in 1:length(toKeep)) {

# 	HuPlex[u] = length(intersect(goHu.bc[[u]], goPlex.bc[[u]]))
# }

# # between A and C
# toKeep = intersect(names(goGrid), names(goPlex))

# goGrid.ac = goGrid[toKeep]
# goPlex.ac = goPlex[toKeep]

# for(u in 1:length(toKeep)) {

# 	GridPlex[u] = length(intersect(goPlex.ac[[u]], goGrid.ac[[u]]))
# }


# sumIntersect = sum(intersection)

# # sumOnlyGrid = sum(onlyGrid)
# # sumOnlyHu = sum(onlyHu)
# # sumOnlyPlex = sum(onlyPlex)

# sumGridHu = sum(GridHu)
# sumHuPlex = sum(HuPlex)
# sumGridPlex = sum(GridPlex)

# draw.triple.venn(sum(lengths(goGrid)), 
# 	sum(lengths(goHu)),
# 	sum(lengths(goPlex)),
# 	sumGridHu,
# 	sumHuPlex,
# 	sumGridPlex,
# 	sumIntersect,
# 	category = c("Original GMT, BioGRID", "Original GMT, huMAP", "Original GMT, BioPlex"),
# 	scaled = FALSE,
# 	fill = c("gold1", "lightskyblue", "lawngreen"),
# 	lty = "blank",
# 	cex = rep(3, 7),
# 	cat.pos = c(-15, 15, 180),
# 	cat.dist = c(0.04, 0.04, 0.03),
# 	cat.cex = c(2, 2, 2)
# 	)


require(data.table)

oG = setkey(data.table(oldGrid))
oH = setkey(data.table(oldHu))
oG.oH = length(na.omit(oH[oG, which=TRUE]))

oH = setkey(data.table(oldHu))
oP = setkey(data.table(oldPlex))
oH.oP = length(na.omit(oP[oH, which=TRUE]))

oG = setkey(data.table(oldGrid))
oP = setkey(data.table(oldPlex))
oG.oP = length(na.omit(oP[oG, which=TRUE]))

int1 = na.omit(oP[oH, which=TRUE])

int2 = length(na.omit(oG[oP[int1, ], which = TRUE]))

draw.triple.venn(nrow(oldGrid), 
	nrow(oldHu),
	nrow(oldPlex),
	oG.oH,
	oH.oP,
	oG.oP,
	int2,
	category = c("Original PPI, BioGRID", "Original PPI, huMAP", "Original PPI, BioPlex"),
	scaled = FALSE,
	fill = c("gold1", "lightskyblue", "lawngreen"),
	lty = "blank",
	cex = rep(3, 7),
	cat.pos = c(-15, 15, 180),
	cat.dist = c(0.04, 0.04, 0.03),
	cat.cex = c(2, 2, 2)
	)