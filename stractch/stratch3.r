data <- read.csv("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/stractch/result/discovery_casecase.csv",stringsAsFactors = F)


try <- data[,2]
try2 <- strsplit(try,"\\(")
try2[[1]][1]
try2[[1]][2]
try3 <- strsplit(try2[[1]][2],"-")
as.numeric(try3[[1]][1])
gsub(")","",try3[[1]][2])
try3 <- strsplit(try3[[1]],"")


data[,2] <- process(data[,2])
data[,4] <- process(data[,4])
data[,6] <- process(data[,6])
data[,8] <- process(data[,8])

process <- function(x){
  places <- 2
  temp.x <- strsplit(x,"\\(")
  n <- length(x)
  result.low <- rep(0,n)
  result <- rep(0,n)
  result.high <- rep(0,n)
  final.result <- rep("c",n)
  for(i in 1:n){
    result[i] <- round(as.numeric(temp.x[[i]][1]),places)
    temp2 <- strsplit(temp.x[[i]][2],"-")
    result.low[i] <- round(as.numeric(temp2[[1]][1]),places)
    result.high[i] <- round(as.numeric(gsub(")","",temp2[[1]][2])),places)
    final.result[i] <- paste0(result[i],"(",result.low[i],"-",result.high[i],")")
    }
  
  return(final.result)
}

write.csv(data,"/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/stractch/result/discovery_casecase_2des.csv")
