#install_github("andrewhaoyu/bcutility",args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
library(data.table)
library(bcutility)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
load("./discovery_SNP/conditional_analysis_check_LD/potential_discovery_snp.Rdata")


discovery_snp_infor <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/SNP_infor_conditional_check.csv",header=T)
discovery_snp_infor <- discovery_snp_infor[1:20,]








library(bc2)
data1 <- as.data.frame(fread("./data/icogs_overall.csv",header=T))
discovery.snp.icog <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis_check_LD/discovery_icog_data.csv"))
country1 <- as.factor(data1[,3])

snp.name <- colnames(discovery.snp.icog)[i1]
########find the snp information
idx <- which(discovery_snp$SNP.ICOGS==snp.name)
snp.chr <- discovery_snp$CHR[idx]
snp.pos <- discovery_snp$position[idx]
########find the corresponding known snp name
idx.infor <- which(discovery_snp_infor$CHR==snp.chr&
                     discovery_snp_infor$Pos==snp.pos)
if(idx.infor==15|idx.infor==19){
  known.snp <- discovery_snp_infor[idx.infor,c(4,5)]
  
}else{
  known.snp <- discovery_snp_infor[idx.infor,4]  
}
idx.known <- which(colnames(data1)%in%known.snp)

x.covar1 <- cbind(data1[,c(5:14)],data1[,idx.known],country1)
gene1 <- discovery.snp.icog[,i1]


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


result <- list(p.value.two.stage.model=p.value.two.stage.model
)

save(result,file=paste0("./discovery_SNP/conditional_analysis_novel/result/novel_conditional_result",i1,".Rdata"))


