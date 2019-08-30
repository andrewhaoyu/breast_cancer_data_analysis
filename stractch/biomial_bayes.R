#n total nodes observed
#x nodes are negative
BionomialBayes <- function(n,x){
  return((x+1)/(n+2))
}
BionomialBayes(0,0)
BionomialBayes(1,1)
BionomialBayes(10,10)
BionomialBayes(10,5)
BionomialBayes(10,1)
