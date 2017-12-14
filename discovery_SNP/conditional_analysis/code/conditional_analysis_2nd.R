rm(list=ls())
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])


source("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/code/conditional_additive_model_update.R")



load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog.clean.sub",i1,".Rdata"))
load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.onco.clean.sub",i1,".Rdata"))

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
library(bigmemory)
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/support.matrix.Rdata")
z.standard <- support.matrix[[1]]
z.additive.design <- support.matrix[[2]]
M <- support.matrix[[3]]
number.of.tumor <- support.matrix[[4]]
z.design.score.baseline <- support.matrix[[5]]
z.design.score.casecase <- support.matrix[[6]]
z.design.score.baseline.ER <- support.matrix[[7]]
z.design.score.casecase.ER <- support.matrix[[8]]



n.condition <- nrow(all.conditional.snps)

start.end <- startend(n.condition,3000,i1)
start <- start.end[1]
end <- start.end[2]

data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
age1 <- data1[,204]
idx.complete1 <- which(age1!=888)
age1 <- age1[idx.complete1]

y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")

known.all.mis1 <- data1[,c(27:203)]
###fake the onco array only snp



x.covar.mis1 <- data1[,c(5:14)]
snp.name.all.icog <- conditional.snp.list.icog.clean.sub[[1]]
x.test.all.mis1 <- conditional.snp.list.icog.clean.sub[[2]]
#x.test.all.mis1 <- conditional.snp.list.icog.clean.sub[[2]]


data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
age2 <- data2[,204]
idx.complete2 <- which(age2!=888)
age2 <- age2[idx.complete2]
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
known.all.mis2 <- data2[,c(27:203,205)]
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")


snp.name.all.onco <- conditional.snp.list.onco.clean.sub[[1]]
x.test.all.mis2 <- conditional.snp.list.onco.clean.sub[[2]]
x.covar.mis2 <- data2[,c(5:14)]



icog.julie <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_icog.csv")
icog.julie <- icog.julie[,-1]
discovery.snp.icog <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T)
onco.julie <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_onco.csv")
onco.julie <- onco.julie[,-1]
discovery.snp.onco <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv")
x.discovery.mis1 <- as.data.frame(cbind(icog.julie,discovery.snp.icog))
x.discovery.mis2 <- as.data.frame(cbind(onco.julie,discovery.snp.onco))
###two snps onco array only
sudo.icog.na <- rep(NA,nrow(data1))

known.all.mis1 <- cbind(known.all.mis1,sudo.icog.na,x.discovery.mis1)
known.all.mis2 <- cbind(known.all.mis2,x.discovery.mis2)




y.pheno.mis1 <- y.pheno.mis1[idx.complete1,]
x.covar.mis1 <- x.covar.mis1[idx.complete1,]
x.covar.mis1 <- cbind(x.covar.mis1,age1)

known.all.mis1 <- known.all.mis1[idx.complete1,]
y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]
x.covar.mis2 <- x.covar.mis2[idx.complete2,]
x.covar.mis2 <- cbind(x.covar.mis2,age2)

known.all.mis2 <- known.all.mis2[idx.complete2,]

known.flag.all <- all.conditional.snps$known.flag



p.value.all <- rep(0,end-start+1)

known.flag.last <- 999
known.flag.new <- 999



load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.first.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.first.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.first")

icog.first.snpvalue <- icog.first[[2]][idx.complete1,]
onco.first.snpvalue <- onco.first[[2]][idx.complete2,]

fine_mapping <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/fine_mapping_annotated_clean.csv",header= T,
                         stringsAsFactors = F)


region.all <- fine_mapping$region.idx
new.region <- c(179:207)
region.all <- c(region.all,new.region)



n.first <- nrow(conditional.results.first)
first.known.flag <- conditional.results.first$known.flag


for(i2 in 1:(end-start+1)){
  print(i2)
  snp.name.icog <- snp.name.all.icog[i2]
  snp.icog <- x.test.all.mis1[,i2]
  snp.onco <- x.test.all.mis2[,i2]
  snp.name.onco <- snp.name.all.onco[i2]
  snp.icog <- snp.icog[idx.complete1]
  snp.onco <- snp.onco[idx.complete2]
  known.flag <- known.flag.all[start+i2-1]
  
  if(known.flag%in%first.known.flag){
    known.flag.new<- known.flag
    if(known.flag.new!=known.flag.last){
      
      load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/score.test.support.icog.2nd",known.flag,".Rdata"))
load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/score.test.support.onco.2nd",known.flag,".Rdata"))
      known.flag.last <- known.flag.new
    }
    
    idx.condition <- which(first.known.flag==known.flag)
    
    conditional.snps.icog <- icog.first.snpvalue[,idx.condition]
    conditional.snps.onco <- onco.first.snpvalue[,idx.condition]
    
    p.value.all[i2] <- condition_additive_model_update(y.pheno.mis1,
                                                x.covar.mis1,
                                                snp.name.icog,
                                                snp.icog,
                                                y.pheno.mis2,
                                                x.covar.mis2,
                                                snp.name.onco,
                                                snp.onco,
                                                known.flag,
                                                known.all.mis1,
                                                known.all.mis2,
                                                z.standard,
                                                z.additive.design,
                                                M,
                                                number.of.tumor,
                                                z.design.score.baseline,
                                                z.design.score.casecase,
                                                z.design.score.baseline.ER,
                                                z.design.score.casecase.ER,
                                                score.test.support.icog = score.test.support.icog,
                                                score.test.support.onco = score.test.support.onco,
                                                conditional.snps.icog=conditional.snps.icog,
                                                conditional.snps.onco=conditional.snps.onco,
                                                region.all = region.all)
    
  }else{
    p.value.all[i2] <- 1
  }
  
  
  
}


save(p.value.all,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/psub.2nd",i1,".Rdata"))

