library(tidyverse)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
source("./discovery_SNP/conditional_analysis/code/all_additive_support_fun.R")

library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/support.matrix.Rdata")
z.standard <- support.matrix[[1]]
z.additive.design <- support.matrix[[2]]
M <- support.matrix[[3]]
number.of.tumor <- support.matrix[[4]]
z.design.score.baseline <- support.matrix[[5]]
z.design.score.casecase <- support.matrix[[6]]
z.design.score.baseline.ER <- support.matrix[[7]]
z.design.score.casecase.ER <- support.matrix[[8]]





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


load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.first.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.first.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.first")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.2nd.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.2nd.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.2nd")



load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.3rd.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.3rd.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.3rd")


load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.4th.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.4th.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.4th")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/icog.5th.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/onco.5th.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.results.5th")



icog.first.snpvalue <- icog.first[[2]][idx.complete1,]
onco.first.snpvalue <- onco.first[[2]][idx.complete2,]

icog.2nd.snpvalue <- icog.2nd[[2]][idx.complete1,]
onco.2nd.snpvalue <- onco.2nd[[2]][idx.complete2,]

icog.3rd.snpvalue <- icog.3rd[[2]][idx.complete1,]
onco.3rd.snpvalue <- onco.3rd[[2]][idx.complete2,]

icog.4th.snpvalue <- icog.4th[[2]][idx.complete1,]
onco.4th.snpvalue <- onco.4th[[2]][idx.complete2,]

icog.5th.snpvalue <- icog.5th[[2]][idx.complete1]
onco.5th.snpvalue <- onco.5th[[2]][idx.complete2]





n.first <- nrow(conditional.results.first)
first.known.flag <- conditional.results.first$known.flag
n.2nd <- nrow(conditional.results.2nd)
known.flag.2nd <- conditional.results.2nd$known.flag

n.3rd <- nrow(conditional.results.3rd)
known.flag.3rd <- conditional.results.3rd$known.flag

n.4th <- nrow(conditional.results.4th)
known.flag.4th <- conditional.results.4th$known.flag

n.5th <- nrow(conditional.results.5th)
known.flag.5th <- conditional.results.5th$known.flag


fine_mapping <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/fine_mapping_annotated_clean.csv",header= T,
                         stringsAsFactors = F)


region.all <- fine_mapping$region.idx
new.region <- c(179:207)
region.all <- c(region.all,new.region)

fine_mapping_dis <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snps_annotated_clean.csv",header= T,
                             stringsAsFactors = F)
sudo.na <- rep(NA,29)
sudo.na.frame <- as.data.frame(matrix(NA,29,9))
fine_mapping_dis_new <- cbind(fine_mapping_dis[,1],sudo.na,fine_mapping_dis[,3],fine_mapping_dis[,2],sudo.na.frame)
colnames(fine_mapping_dis_new) <- colnames(fine_mapping)
fine_mapping <- rbind(fine_mapping,fine_mapping_dis_new)




all.condition.results <- list(conditional.results.first,
                              conditional.results.2nd,
                              conditional.results.3rd,
                              conditional.results.4th,
                              conditional.results.5th)
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata")
##############generate conditional snp list
##############novel snp 189
##############novel snp 190
##############novel snp 194
##############novel snp 195
i1 = 189
idx.first <- which(i1==first.known.flag)
first.snp.value.icog <- icog.first.snpvalue[,idx.first]
first.snp.value.onco <- onco.first.snpvalue[,idx.first]
idx.known <- which(region.all==region.all[i1])  
known.snp.value.icog <- as.matrix(known.all.mis1[,idx.known])
snp.icog <- known.snp.value.icog[,1]
known.snp.value.onco <- as.matrix(known.all.mis2[,idx.known])
snp.onco <-   known.snp.value.onco[,1]
idx.2nd <- which(i1==known.flag.2nd)
snp.value.icog.2nd <- icog.2nd.snpvalue[,idx.2nd]
snp.value.onco.2nd <- onco.2nd.snpvalue[,idx.2nd]
idx.3rd <- which(i1==known.flag.3rd)
snp.value.icog.3rd <- icog.3rd.snpvalue[,idx.3rd]
snp.value.onco.3rd <- onco.3rd.snpvalue[,idx.3rd]
idx.4th <- which(i1==known.flag.4th)
snp.value.icog.4th <- icog.4th.snpvalue[,idx.4th]
snp.value.onco.4th <- onco.4th.snpvalue[,idx.4th]
idx.5th <- which(i1==known.flag.5th)
snp.value.icog.5th <- icog.5th.snpvalue
snp.value.onco.5th <- onco.5th.snpvalue
all.idx <- list(idx.first,idx.2nd,idx.3rd,idx.4th,idx.5th)
idx.near.known.icog <- which(colnames(data1)=="rs11249433")
near.known.icog <- data1[idx.complete1,idx.near.known.icog]
idx.near.known.onco <- which(colnames(data2)=="rs11249433")
near.known.onco <- data2[idx.complete2,idx.near.known.onco]


x.all.mis.icog.1 <- cbind(near.known.icog,known.snp.value.icog,first.snp.value.icog,snp.value.icog.2nd,snp.value.icog.3rd,snp.value.icog.4th,snp.value.icog.5th)
x.all.mis.onco.1 <- cbind(near.known.onco,known.snp.value.onco,first.snp.value.onco,snp.value.onco.2nd,snp.value.onco.3rd,snp.value.onco.4th,snp.value.onco.5th)

i1 = 190


idx.first <- which(i1==first.known.flag)
first.snp.value.icog <- icog.first.snpvalue[,idx.first]
first.snp.value.onco <- onco.first.snpvalue[,idx.first]
idx.known <- which(region.all==region.all[i1])  
known.snp.value.icog <- as.matrix(known.all.mis1[,idx.known])
snp.icog <- known.snp.value.icog[,1]
known.snp.value.onco <- as.matrix(known.all.mis2[,idx.known])
snp.onco <-   known.snp.value.onco[,1]
idx.2nd <- which(i1==known.flag.2nd)
snp.value.icog.2nd <- icog.2nd.snpvalue[,idx.2nd]
snp.value.onco.2nd <- onco.2nd.snpvalue[,idx.2nd]
idx.3rd <- which(i1==known.flag.3rd)
snp.value.icog.3rd <- icog.3rd.snpvalue[,idx.3rd]
snp.value.onco.3rd <- onco.3rd.snpvalue[,idx.3rd]
idx.near.known.icog <- which(colnames(data1)=="rs10941679")
near.known.icog <- data1[idx.complete1,idx.near.known.icog]
idx.near.known.onco <- which(colnames(data2)=="rs10941679")
near.known.onco <- data2[idx.complete2,idx.near.known.onco]


x.all.mis.icog.2 <- cbind(near.known.icog,known.snp.value.icog,first.snp.value.icog,snp.value.icog.2nd,snp.value.icog.3rd)
x.all.mis.onco.2 <- cbind(near.known.onco,known.snp.value.onco,first.snp.value.onco,snp.value.onco.2nd,snp.value.onco.3rd)

i1 = 194


idx.first <- which(i1==first.known.flag)
first.snp.value.icog <- icog.first.snpvalue[,idx.first]
first.snp.value.onco <- onco.first.snpvalue[,idx.first]
idx.known <- which(region.all==region.all[i1])  
known.snp.value.icog <- as.matrix(known.all.mis1[,idx.known])
snp.icog <- known.snp.value.icog[,1]
known.snp.value.onco <- as.matrix(known.all.mis2[,idx.known])
snp.onco <-   known.snp.value.onco[,1]
idx.2nd <- which(i1==known.flag.2nd)
snp.value.icog.2nd <- icog.2nd.snpvalue[,idx.2nd]
snp.value.onco.2nd <- onco.2nd.snpvalue[,idx.2nd]
idx.3rd <- which(i1==known.flag.3rd)
snp.value.icog.3rd <- icog.3rd.snpvalue[,idx.3rd]
snp.value.onco.3rd <- onco.3rd.snpvalue[,idx.3rd]
idx.near.known.icog <- which(colnames(data1)=="rs10941679")
near.known.icog <- data1[idx.complete1,idx.near.known.icog]
idx.near.known.onco <- which(colnames(data2)=="rs10941679")
near.known.onco <- data2[idx.complete2,idx.near.known.onco]


x.all.mis.icog.3 <- cbind(near.known.icog,known.snp.value.icog,first.snp.value.icog,snp.value.icog.2nd,snp.value.icog.3rd)
x.all.mis.onco.3 <- cbind(near.known.onco,known.snp.value.onco,first.snp.value.onco,snp.value.onco.2nd,snp.value.onco.3rd)

i1 = 195
idx.first <- which(i1==first.known.flag)
first.snp.value.icog <- icog.first.snpvalue[,idx.first]
first.snp.value.onco <- onco.first.snpvalue[,idx.first]
idx.known <- which(region.all==region.all[i1])  
known.snp.value.icog <- as.matrix(known.all.mis1[,idx.known])
snp.icog <- known.snp.value.icog[,1]
known.snp.value.onco <- as.matrix(known.all.mis2[,idx.known])
snp.onco <-   known.snp.value.onco[,1]
idx.near.known.icog <- which(colnames(data1)%in%c("rs17879961","rs132390"))
near.known.icog <- data1[idx.complete1,idx.near.known.icog]
idx.near.known.onco <- which(colnames(data2)%in%c("rs17879961","rs132390"))
near.known.onco <- data2[idx.complete2,idx.near.known.onco]

x.all.mis.icog.4 <- cbind(near.known.icog,known.snp.value.icog,first.snp.value.icog)
x.all.mis.onco.4 <- cbind(near.known.onco,known.snp.value.onco,first.snp.value.onco)


conditional.check.data <- list(x.all.mis.icog.1,
                               x.all.mis.onco.1,
                               x.all.mis.icog.2,
                               x.all.mis.onco.2,
                               x.all.mis.icog.3,
                               x.all.mis.onco.3,
                               x.all.mis.icog.4,
                               x.all.mis.onco.4)
save(conditional.check.data,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.check.data.Rdata")