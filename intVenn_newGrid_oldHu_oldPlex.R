# say we have goHu, goPlex, goGrid - the old GMTs
# newGmt signifies the new GMTs

# Part II: trying to find common interactions between new GRID, old Plex, and old huMAP

# Part I: Do common interactions between all 3.

# gives us all the common protein sets between all three (new, old, old)
# toKeep = intersect(intersect(names(newGmtGrid), names(goHu)), names(goPlex))

# # gives us three GMTs in list format with the same protein set names in the same order
# newGmtGrid.abc = newGmtGrid[toKeep]
# goHu.abc = goHu[toKeep]
# goPlex.abc = goPlex[toKeep]

# intersection = c() # intersection of all three

# for (u in 1:length(toKeep)) {

# 	intersection[u] = length(intersect(intersect(newGmtGrid.abc[[u]], goHu.abc[[u]]), goPlex.abc[[u]]))

# }

# # Part II: Do common interactions between pairs of 2, including the interactions between all 3.

# GridHu = c() # both new BioGRID and old huMAP (including all 3)
# HuPlex = c() # both old BioPlex and old huMAP
# GridPlex = c() # both new BioGRID and old BioPlex

# # between A and B
# toKeep = intersect(names(newGmtGrid), names(goHu))

# newGmtGrid.ab = newGmtGrid[toKeep]
# goHu.ab = goHu[toKeep]

# for (u in 1:length(toKeep)) {

# 	GridHu[u] = length(intersect(newGmtGrid.ab[[u]], goHu.ab[[u]]))
# }

# # between B and C
# toKeep = intersect(names(goHu), names(goPlex))

# goHu.bc = goHu[toKeep]
# goPlex.bc = goPlex[toKeep]

# for (u in 1:length(toKeep)) {

# 	HuPlex[u] = length(intersect(goHu.bc[[u]], goPlex.bc[[u]]))
# }

# # between A and C
# toKeep = intersect(names(newGmtGrid), names(goPlex))

# newGmtGrid.ac = newGmtGrid[toKeep]
# goPlex.ac = goPlex[toKeep]

# for(u in 1:length(toKeep)) {

# 	GridPlex[u] = length(intersect(goPlex.ac[[u]], newGmtGrid.ac[[u]]))
# }


# sumIntersect = sum(intersection)

# sumGridHu = sum(GridHu)
# sumHuPlex = sum(HuPlex)
# sumGridPlex = sum(GridPlex)

# draw.triple.venn(sum(lengths(newGmtGrid)), 
# 	             sum(lengths(goHu)),
# 	             sum(lengths(goPlex)),
# 	             sumGridHu,
# 	             sumHuPlex,
# 	             sumGridPlex,
# 	             sumIntersect,
# 	             category = c("New GMT, BioGRID", "Original GMT, huMAP", "Original GMT, BioPlex"),
# 			 scaled = FALSE,
# 			 fill = c("goldenrod3", "lightskyblue", "lawngreen"),
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
oP = setkey(data.table(oldPlex))
oH.oP = length(na.omit(oP[oH, which=TRUE]))

nG = setkey(data.table(newGrid))
oP = setkey(data.table(oldPlex))
nG.oP = length(na.omit(oP[nG, which=TRUE]))

inter = na.omit(oP[oH, which=TRUE])

intersection = length(na.omit(nG[oP[inter, ], which = TRUE]))

draw.triple.venn(nrow(newGrid), 
	nrow(oldHu),
	nrow(oldPlex),
	nG.oH,
	oH.oP,
	nG.oP,
	intersection,
	category = c("New PPI, BioGRID", "Original PPI, huMAP", "Original PPI, BioPlex"),
	scaled = FALSE,
	fill = c("darkgoldenrod1", "lightskyblue", "lawngreen"),
	lty = "blank",
	cex = rep(3, 7),
	cat.pos = c(-15, 15, 180),
	cat.dist = c(0.04, 0.04, 0.03),
	cat.cex = c(2, 2, 2)
	)