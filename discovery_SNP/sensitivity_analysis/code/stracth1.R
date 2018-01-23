y <- y.pheno.complete
tumor.number <- ncol(y)-1
y.case.control <- y[,1]
y.tumor <- y[,2:(tumor.number+1)]
freq.subtypes <- GenerateFreqTable(y.pheno.complete)

freq.subtypes[freq.subtypes[,5]>10,]






idx.clean <- NULL
rowmatch <- function(A,B) { 
  # Rows in A that match the rows in B 
  f <- function(...) paste(..., sep=":") 
  if(!is.matrix(B)) B <- matrix(B, 1, length(B)) 
  a <- do.call("f", as.data.frame(A)) 
  b <- do.call("f", as.data.frame(B)) 
  # a <- apply(A,1,paste0,collapse=":")
  which(a %in% b)
} 

idx.temp <- rowmatch(y[,2:ncol(y)],subtype.result)


