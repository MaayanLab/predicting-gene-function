get_ttest <- function(data) {
  
  # split up the violin data matrices
  gtex1 = match("GTEX", data[, 3])
  hum1 = match("Human", data[, 3])
  ccle1 = match("CCLE", data[, 3])
  
  # for Mouse to GTEX
  statG = t.test(data[1:(gtex1 - 1), 2], data[gtex1:(hum1 - 1), 2])$statistic
  pValG = t.test(data[1:(gtex1 - 1), 2], data[gtex1:(hum1 - 1), 2])$p.value
  dAvG = t.test(data[1:(gtex1 - 1), 2], data[gtex1:(hum1 - 1), 2])$estimate[1] - 
         t.test(data[1:(gtex1 - 1), 2], data[gtex1:(hum1 - 1), 2])$estimate[2]
  
  # for Mouse to Human
  statH = t.test(data[1:(gtex1 - 1), 2], data[hum1:(ccle1 - 1), 2])$statistic
  pValH = t.test(data[1:(gtex1 - 1), 2], data[hum1:(ccle1 - 1), 2])$p.value
  dAvH = t.test(data[1:(gtex1 - 1), 2], data[hum1:(ccle1 - 1), 2])$estimate[1] - 
         t.test(data[1:(gtex1 - 1), 2], data[hum1:(ccle1 - 1), 2])$estimate[2]
  
  # for Mouse to CCLE
  statC = t.test(data[1:(gtex1 - 1), 2], data[ccle1:nrow(data), 2])$statistic
  pValC = t.test(data[1:(gtex1 - 1), 2], data[ccle1:nrow(data), 2])$p.value
  dAvC = t.test(data[1:(gtex1 - 1), 2], data[ccle1:nrow(data), 2])$estimate[1] - 
         t.test(data[1:(gtex1 - 1), 2], data[ccle1:nrow(data), 2])$estimate[2]
  
  results = matrix(c(statG, pValG, dAvG, statH, pValH, dAvH, statC, pValC, dAvC), nrow = 3, ncol = 3,
                   dimnames = list(c("Statistic", "p-value", paste0(expression(delta), " Mean")), 
                                   c("Mouse to GTEX", "Mouse to Human", "Mouse to CCLE")))
  
  return(results)
  
}

for (a in 1:length(listVd)) {
  
  write.csv(get_ttest(listVd[[a]]), file= paste0(names(listVd)[a], "ttest.csv"))
  
}