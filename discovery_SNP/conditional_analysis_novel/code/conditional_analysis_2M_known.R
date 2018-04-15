#install_github("andrewhaoyu/bcutility",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
library(data.table)
library(bcutility)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
discovery_snp <- as.data.frame(fread("./data/discovery_snps_annotated_clean.csv",header=T))
fine_mapping <- as.data.frame(fread("./data/fine_mapping_annotated_clean.csv"))
library(bc2)
dis <- 2*10^6
check.data <- NULL
check.id <- NULL
for(i in 11:28){
  pos <- discovery_snp$position[i]
  CHR <- discovery_snp$CHR[i]
  for(j in 1:178){
    pos.known <- fine_mapping$position[j]
    CHR.known <- fine_mapping$CHR[j]
    if(CHR==CHR.known&pos>=(pos.known-dis)&pos<=(pos.known+dis)){
      print(c(i-10,j))
      check.id <- rbind(check.id,c(i,j))
      temp1 <- discovery_snp[i,c(1,3,2)]
      colnames(temp1) <- c("SNP","CHR","Position")
      temp2 <- fine_mapping[j,c(1,3,4)]
      colnames(temp2) <- c("SNP","CHR","Position")
      result <- rbind(temp1,temp2)
      check.data <- rbind(check.data,result)      
    }
    
  }
}

#write.csv(check.data,file="./data/check_SNPs.csv",row.names = F,col.names = T)

data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
discovery.snp.icog <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T))
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
# Grade1.fake <- data1$Grade1
# Grade1.fake[data1$Grade1==2|data1$Grade1==3] <- 1
# Grade1.fake[data1$Grade1==1] <- 0
#y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1,Grade1.fake)
# y.pheno.mis1 <- cbind(data1$Behaviour1,data1$PR_status1,data1$ER_status1,data1$HER2_status1)

#x.test.all.mis1 <- data1[,c(27:203)]
###pc1-10 and age
#i1 <- 4

idx.check <- which(check.data[,1]==check.data[2*i1-1,1])
if(length(idx.check)==1){
  idx.known <- which(colnames(data1)==check.data[idx.check+1,1])  
}else{
  idx.known <- NULL
  for(i in 1:length(idx.check)){
    idx.temp <- which(colnames(data1)==check.data[idx.check[i]+1,1])  
    idx.known <- c(idx.known,idx.temp)
  }
}


x.covar1 <- cbind(data1[,c(5:14)],data1[,idx.known])

gene1 <- discovery.snp.icog[,check.id[i1,1]-10]
age <- data1[,204]
x.covar1 <- cbind(x.covar1,age)
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar1 <- x.covar1[idx.complete,]
gene1 <- gene1[idx.complete]


data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv"))

data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

# x.test.all.mis2 <- data2[,c(27:203)]
if(length(idx.check)==1){
  idx.known <- which(colnames(data2)==check.data[idx.check+1,1])  
}else{
  idx.known <- NULL
  for(i in 1:length(idx.check)){
    idx.temp <- which(colnames(data2)==check.data[idx.check[i]+1,1])  
    idx.known <- c(idx.known,idx.temp)
  }
}

x.covar2 <- cbind(data2[,c(5:14)],data2[,idx.known])
age <- data2[,204]
x.covar2 <- cbind(x.covar2,age)
idx.complete <- which(age!=888)
gene2 <- discovery.snp.onco[,check.id[i1,1]-10]
y.pheno.mis2 <- y.pheno.mis2[idx.complete,]
x.covar2 <- x.covar2[idx.complete,]
gene2 <- gene2[idx.complete]
#gene.known2 <- x.covar2[,11]
#idx.control <- which(y.pheno.mis2[,1]==0)
#cor(gene2[idx.control],gene.known2[idx.control])^2
z.standard <- GenerateZstandard(y.pheno.mis1)
z.random.support <- cbind(1,z.standard[,1])
z.random.test <- z.standard[,2:4]
p.value.two.stage.model <- two_data_two_stage_random(y.pheno.mis1,
                                                     gene1,
                                                     x.covar1,
                                                     y.pheno.mis2,
                                                     gene2,
                                                     x.covar2)

data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
discovery.snp.icog <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T))
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")

country1 <- as.factor(data1[,3])
x.covar1 <- cbind(data1[,c(5:14)],data1[,idx.known],country1)
gene1 <- discovery.snp.icog[,check.id[i1,1]-10]

idx.check <- which(check.data[,1]==check.data[2*i1-1,1])
if(length(idx.check)==1){
  idx.known <- which(colnames(data1)==check.data[idx.check+1,1])  
}else{
  idx.known <- NULL
  for(i in 1:length(idx.check)){
    idx.temp <- which(colnames(data1)==check.data[idx.check[i]+1,1])  
    idx.known <- c(idx.known,idx.temp)
  }
}

data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv"))

data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

# x.test.all.mis2 <- data2[,c(27:203)]
if(length(idx.check)==1){
  idx.known <- which(colnames(data2)==check.data[idx.check+1,1])  
}else{
  idx.known <- NULL
  for(i in 1:length(idx.check)){
    idx.temp <- which(colnames(data2)==check.data[idx.check[i]+1,1])  
    idx.known <- c(idx.known,idx.temp)
  }
}
country2 <- as.factor(data2[,4])
x.covar2 <- cbind(data2[,c(5:14)],data2[,idx.known],country2)
gene2 <- discovery.snp.onco[,check.id[i1,1]-10]
#age <- data2[,204]

p.value.standard <- two_data_standard_anlysis(y.pheno.mis1,
                                              gene1,
                                              x.covar1,
                                              y.pheno.mis2,
                                              gene2,
                                              x.covar2)


result <- list(p.value.two.stage.model=p.value.two.stage.model,
               p.value.standard= p.value.standard)

save(result,file=paste0("./discovery_SNP/conditional_analysis_novel/novel_conditional_reuslt",i1,".Rdata"))


