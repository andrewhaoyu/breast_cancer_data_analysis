#install_github("andrewhaoyu/bcutility",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
library(data.table)
library(bcutility)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
load("./discovery_SNP/conditional_analysis_check_LD/potential_discovery_snp.Rdata")


discovery_snp_infor <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/SNP_infor_conditional_check.csv",header=T,stringsAsFactors = F)
discovery_snp_infor <- discovery_snp_infor[1:20,]








library(bc2)
data1 <- as.data.frame(fread("./data/icogs_overall.csv",header=T))
idx.case1 <- which(data1$Behaviour1==2|data1$Behaviour1==888)
discovery.snp.icog <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis_check_LD/discovery_icog_data_overall.csv"))
colnames(discovery.snp.icog)[5] <- "5:45333860"
country1 <- as.factor(data1[,3])
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
y.pheno.mis1[idx.case1,1] <- 1
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis1) = c("Behaviour","ER",
                           "PR","HER2","Grade")

snp.name <- colnames(discovery.snp.icog)[i1]
#################find the known snp information
temp <- strsplit(snp.name,":")
pos <- as.numeric(temp[[1]][2])
idx.dis.infor <- which(discovery_snp_infor$Pos==pos)
if(discovery_snp_infor[idx.dis.infor,5]==""){
  near.known.snp.name <- discovery_snp_infor[idx.dis.infor,4]  
}else{
  near.known.snp.name <- discovery_snp_infor[idx.dis.infor,c(4,5)]  
}
#idx.known1 <- which(colnames(data1)%in%near.known.snp.name)


#idx.known <- which(colnames(data1)%in%known.snp)
x.covar1 <- cbind(data1[,c(6:15)],country1)
#x.covar1 <- cbind(data1[,c(6:15)],data1[,idx.known1],country1)
gene1 <- discovery.snp.icog[,i1]

%>% %>% 
data2 <- fread("./data/oncoarray_overall.csv",header=T)
discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis_check_LD/discovery_onco_data_overall.csv"))

data2 <- as.data.frame(data2)
idx.case2 <- which(data2$Behaviour1==2|data2$Behaviour1==888)
#data2$Behaviour1[idx.case2] <- 1
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
y.pheno.mis2[idx.case2,1] <- 1
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")

# x.test.all.mis2 <- data2[,c(27:203)]
#idx.known2 <- which(colnames(data2)%in%near.known.snp.name)
country2 <- as.factor(data2[,4])
x.covar2 <- cbind(data2[,c(5:14)],country2)
#x.covar2 <- cbind(data2[,c(5:14)],data2[,idx.known2],country2)
gene2 <- discovery.snp.onco[,i1]
#age <- data2[,204]


p.value.standard <- two_data_standard_anlysis(y.pheno.mis1,
                                              gene1,
                                              x.covar1,
                                              y.pheno.mis2,
                                              gene2,
                                              x.covar2)

result <- list(p.value.standard,snp.name
)

save(result,file=paste0("./discovery_SNP/conditional_analysis_novel/result/standard_conditional_result",i1,".Rdata"))


