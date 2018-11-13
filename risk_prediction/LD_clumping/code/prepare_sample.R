library(data.table)
subject.file <- "/spin1/users/zhangh24/test/onco_order.txt"

onco.order <- read.table(subject.file)
n <- nrow(onco.order)
ID <- matrix("c",n,1)
for(i in 1:n){
  ID[i] <- paste0("sample_",onco.order[i,1])
}
missing <- matrix(0,n,1)
case <- matrix(rbinom(n,1,0.5),n,1)
cov_1 <- matrix(rnorm(n),n,1)
#onco.order <- matrix(paste0("sample",onco.order),n,1)
ID <- rbind(0,ID)
missing <- rbind(0,missing)
case <- rbind("B",case)
cov_1 <- rbind("C",cov_1)
#onco.order <- matrix(onco.order,ncol=1)

sample.data <- data.frame(ID,
                          ID,
                          missing,
                          case,
                          cov_1,stringsAsFactors = T)
colnames(sample.data) <- c("ID_1",
                                "ID_2",
                                "missing",
                                "case",
                                "cov_1")

write.table(sample.data,file = "/spin1/users/zhangh24/test/sample.txt",
            row.names = F, quote = F,sep = " ")

