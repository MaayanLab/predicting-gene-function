
# puts gmt text file into list form
toList <- function(txtf) {

gmt = readLines(txtf)

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

# puts two col form into a text file
# arguments are the PPI in two column form (matrix) and filename (including ".txt")
toPpiTxtFile <- function(twoColForm, filename) {

  sink(filename)

  for (d in 1:nrow(twoColForm)) {
    
    cat(twoColForm[d, 1], twoColForm[d, 2], sep = "\t\t\t\t\t")
    cat("\n")
  }

  sink()

}

# takes ppi in set form and turns it into a two column matrix of paired interactions
toTwoCol <- function(setForm) {

  lst.names = c(1:sum(lengths(setForm)))
  twoCol = vector("list", length(lst.names))
  names(twoCol) = lst.names

  ct = 1

  for (i in 1:length(setForm)) {

    for (j in 1:length(setForm[[i]])) {

      if (names(setForm)[i] != (setForm[[i]])[j]) {

        twoCol[[ct]] = c(names(setForm)[i], (setForm[[i]])[j])
        ct = ct + 1

      }

    }
  }

  # make sure everything is in same order and remove duplicates
  twoCol = unique(lapply(twoCol, sort))

  lth = length(twoCol)

  twoCol = matrix(unlist(twoCol), nrow = lth, ncol = 2, byrow = T)
  return(twoCol)
}