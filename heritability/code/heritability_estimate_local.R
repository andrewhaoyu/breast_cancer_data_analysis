#-------------------------------------------------------------------
# Goal: estimate heritability for breast cancer overall risk and subtypes risk on local
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------

SigmaEst <- function(data){
  #first column is log odds ratio
  #second column is variance
  #third column is allele frequency
  idx <- which(is.na(data[,1]))
  if(length(idx!=0)){
    data <- data[-idx,]  
  }
  
  beta <- data[,1]
  sigma <- data[,2]
  p <- data[,3]
  
  n <- length(p)
  result <- sum(2*p*(1-p)*(beta^2-sigma))
  return(result)
}

LD_pruning = function(sig_SNPs,LD2){
  sig_SNPs_temp =sig_SNPs
  filter_result = NULL
  LD2.temp =LD2
  temp.ind = 1
  while(nrow(sig_SNPs_temp)!=0){
    
    idx = which.min(sig_SNPs_temp$p.value)
    filter_result = rbind(filter_result,sig_SNPs_temp[idx,])
    LD2.single = LD2.temp[idx,]
    idx.cut = which( LD2.single>=0.1)
    position.range <- 500*10^3
    filter_result_position = sig_SNPs_temp$position[idx]
    filter_CHR = sig_SNPs_temp$CHR[idx]
    idx.cut2 <- which((sig_SNPs_temp$position>=filter_result_position-position.range)&(sig_SNPs_temp$position<=filter_result_position+position.range)&(sig_SNPs_temp$CHR==filter_CHR ))
    idx.cut <- c(idx.cut,idx.cut2)
    idx.cut <- unique(idx.cut)
    LD2.temp = LD2.temp[-idx.cut,-idx.cut]
    LD2.temp = as.matrix(LD2.temp)
    sig_SNPs_temp = sig_SNPs_temp[-idx.cut,]
    temp.ind = temp.ind+1
  }
  return(filter_result)
}



library(data.table)
setwd('/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis')
#load summary level statistics for subtypes
load("/dcl01/chatterj/data/hzhang1/breast_intrinsic/whole_genome_breast_cancer_results/BCAC_subtypes_result.Rdata")
#load summary level statistics for overall
standard_result <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF_Haoyuupdate.txt"))
#generate the summary level statistics based on iCOGS and ONCOArray
# Twometa <- function(beta1,var1,beta2,var2){
#   var_meta <- 1/(1/var1+1/var2)
#   beta_meta <- (var_meta)*(beta1/var1+
#                                beta2/var2)
#   return(list(beta_meta,var_meta))
# }
# library(dplyr)
# 
# 
# standard_result = standard_result %>% 
#   mutate(BCAC_meta_beta = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[1]],
#          BCAC_meta_var = Twometa(beta.iCOGs,SE.iCOGs^2,beta.Onco,SE.Onco^2)[[2]],
#          Z = BCAC_meta_beta/sqrt(BCAC_meta_var),
#          P_BCAC_meta = 2*pnorm(-abs(Z))
#          )
# write.table(standard_result,file = "./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF_Haoyuupdate.txt",row.names = F,
# col.names=T,quote=F)



##################discovery snp were ordered based on the order they are extracted
discovery_snp <- read.csv("./data/discovery_snp_summary_new.csv",header=T,stringsAsFactors = F)


discovery_snp_paper_order <- read.csv("./data/discovery_snp_paper_order.csv",header=T,stringsAsFactors = F)
chr.pos.paper <- paste0(discovery_snp_paper_order$CHR,":",discovery_snp_paper_order$position)
library(dplyr)

chr.pos <- paste0(discovery_snp$CHR.x,":",discovery_snp$position)

idx.match <- match(chr.pos.paper,
                   chr.pos)
discovery_snp_new <- discovery_snp[idx.match,]
##create subtypes summary level statistics for 32 new discovery SNPs
BCAC_subtypes_result = BCAC_subtypes_result %>% 
  mutate(chr.pos = paste0(CHR,":",position))
discovery_snp_subtypes <- left_join(discovery_snp_new,BCAC_subtypes_result,by="chr.pos")
##create standard summary level statistics for 32 new discovery SNPs
standard_result = standard_result %>% 
  mutate(chr.pos = paste0(chr.Onco,":",Position.Onco))
discovery_snp_standard <- left_join(discovery_snp_new,standard_result,by="chr.pos")
##create subtypes summary level statistics for 178 known SNPs
fine_map <- read.csv("./data/fine_mapping_annotated_clean.csv")
fine_map <- fine_map %>% 
  mutate(chr.pos = paste0(CHR,":",position))
known_snp_subtypes <- left_join(fine_map,BCAC_subtypes_result,
                            by="chr.pos")
#two snps duplicated
known_snp_subtypes <- known_snp_subtypes[-c(62,120),]
##create standard summary level statistics for 178 known SNPs
known_snp_standard <- left_join(fine_map,standard_result,
                            by="chr.pos")
known_snp_standard <- known_snp_standard[-c(62,120),]

#LD pruning for known SNPs
data2 <- as.data.frame(fread("./data/Onco_euro_v10_10232017.csv",header=T))
names2 = colnames(data2)
idxi1 = which(names2=="rs554219")
x.test.all.mis2 <- data2[,c(27:203,idxi1)]

p.value.known <- known_snp_standard$P_BCAC_meta

idx.control <- which(data2$Behaviour1==0)

LD2 <- cor(x.test.all.mis2[idx.control,])^2
known.snp.infor.p <- known_snp_standard %>% 
  select(Best.published.SNP,CHR,position,P_BCAC_meta)
  
colnames(known.snp.infor.p) <- c("SNP","CHR","position",
                                 "p.value")

known.snp.infor.pruned <- LD_pruning(known.snp.infor.p,LD2)

idx.fil <- which(known.snp.infor.p$SNP%in%known.snp.infor.pruned$SNP)
#only keep SNPs after LD pruning
known_snp_subtypes <- known_snp_subtypes[idx.fil,]
known_snp_standard <- known_snp_standard[idx.fil,]

known_snp_standard_cal = known_snp_standard %>% 
  select(BCAC_meta_beta,BCAC_meta_var,EAFcontrols.Onco)

#heritability for known snps
SigmaEst(known_snp_standard_cal)
#heritability for discovery snps
discovery_snp_cal <- discovery_snp_standard %>% 
  select(BCAC_meta_beta,BCAC_meta_var,EAFcontrols.Onco)
SigmaEst(discovery_snp_cal)
SigmaEst(known_snp_standard_cal)+SigmaEst(discovery_snp_cal)

#heritability for known SNPs subtypes
idx.known <- grep("log_or",colnames(known_snp_subtypes))
jdx.known <- grep("var_",colnames(known_snp_subtypes))
zdx.known <- grep("MAF",colnames(known_snp_subtypes))
idx.dis<- grep("log_or",colnames(discovery_snp_subtypes))
jdx.dis <- c(32:36)
zdx.dis <- grep("exp_freq_a1.x",colnames(discovery_snp_subtypes))
result.known <- rep(0,length(idx))
result.discovery <- rep(0,length(jdx))
for(k in 1:length(idx)){
  known_snp_subtypes_cal <- known_snp_subtypes[,c(idx.known[k],jdx.known[k],zdx.known)]
  result.known[k] <-  SigmaEst(known_snp_subtypes_cal)
  discovery_snp_cal <- discovery_snp_subtypes[,
                                              c(idx.dis[k],jdx.dis[k],zdx.dis)]
  result.discovery[k] <-  SigmaEst(discovery_snp_cal)
}

result.known+result.discovery


#heritability for discovery SNPs subtypes





































result_standard <- cbind(as.numeric(discovery_snp_new$bcac_onco2_beta),
                         as.numeric(discovery_snp_new$bcac_onco2_se)^2)
#discovery_snp_new %>% 
#select(bcac_onco_icogs_gwas_beta,
#      bcac_onco_icogs_gwas_se)
###convert the result to log odds ratio
n.snp <- 32
n.subtypes <- 5
dis.result <- result
# for(i in 1:n.snp){
#   for(j in c(2*(1:n.subtypes)-1)){
#     dis.result[i,j] <- log(as.numeric(strsplit(result[i,j],"\\(")[[1]][1]))
#     
#   }
# }
# dis.result[,c(2*(1:n.subtypes))] <-as.matrix(result[,c(2*(1:n.subtypes))])
# dis.result <- as.data.frame(dis.result)
dis.result <- as.data.frame(dis.result)
colnames(dis.result)[c(2*(1:n.subtypes))-1] <- paste0("logodds_",c("Luminial_A","Luminal_B",
                                                                   "Luminal_B_HER2-",
                                                                   "HER2_Enriched",
                                                                   "Triple_Negative"))
colnames(dis.result)[c(2*(1:n.subtypes))] <- paste0("var_",c("Luminial A","Luminal B",
                                                             "Luminal B HER2-",
                                                             "HER2 Enriched",
                                                             "Triple Negative"))
colnames(result_standard) <- c("logoddsd_standard",
                               "var_standard")
dis.result <- cbind(dis.result,result_standard,freq)
#dis.result[,c(11,12)] <- result_standard






##load the results for 178 known SNPs
known.result <- matrix(0,178,13)
colnames(known.result) <- colnames(dis.result)


for(i in 1:178){
  print(i)
  load(paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_origin",i,".Rdata"))
  known.result[i,] <- Outvec(heter.result)
}
#Onco array paper doens't have the 178th known SNP 
temp.result <- known.result[178,c(11:12)]
###load Onco array results
fine_map <- read.csv("./data/fine_mapping_annotated_clean.csv")
fine_map <- fine_map %>% 
  mutate(chr.pos = paste0(CHR,":",position))
fine_map_new <- left_join(fine_map,all_result,by="chrpos")
#one snp duplicated
fine_map_new <- fine_map_new[-c(62,120),]
known.result[,c(11,12)] <- cbind(as.numeric(fine_map_new$bcac_onco2_beta),as.numeric(fine_map_new$bcac_onco2_se)^2)
known.result[178,c(11,12)] <- temp.result

#########load onco array dataset
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$ER_status1,data2$PR_status1,data2$HER2_status1,data2$Grade1)
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","ER",
                           "PR","HER2","Grade")
names2 = colnames(data2)
idxi1 = which(names2=="rs554219")
x.test.all.mis2 <- data2[,c(27:203,idxi1)]
#LD pruning and position pruning
known.snp.infor <- fine_map[,c(1,3,4)]
discovery.snp.infor <- discovery_snp_new[,c(3,4,5)]
colnames(known.snp.infor) <- colnames(discovery.snp.infor) <- 
  c("SNP","chr","position")
all.snp.infor <- rbind(known.snp.infor,discovery.snp.infor)
x.test.all <- cbind(x.test.all.mis2,x.test.all.mis.dis)
colnames(x.test.all) <- all.snp.infor[,1]
p.value.known <- fine_map_new$bcac_onco_icogs_gwas_P1df

idx.control <- which(y.pheno.mis2[,1]==0)

LD2 <- cor(x.test.all.mis2[idx.control,])^2
known.snp.infor.p <- cbind(known.snp.infor,p.value.known)
colnames(known.snp.infor.p) <- c("SNP","CHR","position",
                                 "p.value")
known.snp.infor.pruned <- LD_pruning(known.snp.infor.p,LD2)

idx.fil <- which(known.snp.infor.p$SNP%in%known.snp.infor.pruned$SNP)









LD_pruning = function(sig_SNPs,LD2){
  sig_SNPs_temp =sig_SNPs
  filter_result = NULL
  LD2.temp =LD2
  temp.ind = 1
  while(nrow(sig_SNPs_temp)!=0){
    
    idx = which.min(sig_SNPs_temp$p.value)
    filter_result = rbind(filter_result,sig_SNPs_temp[idx,])
    LD2.single = LD2.temp[idx,]
    idx.cut = which( LD2.single>=0.1)
    position.range <- 500*10^3
    filter_result_position = sig_SNPs_temp$position[idx]
    filter_CHR = sig_SNPs_temp$CHR[idx]
    idx.cut2 <- which((sig_SNPs_temp$position>=filter_result_position-position.range)&(sig_SNPs_temp$position<=filter_result_position+position.range)&(sig_SNPs_temp$CHR==filter_CHR ))
    idx.cut <- c(idx.cut,idx.cut2)
    idx.cut <- unique(idx.cut)
    LD2.temp = LD2.temp[-idx.cut,-idx.cut]
    LD2.temp = as.matrix(LD2.temp)
    sig_SNPs_temp = sig_SNPs_temp[-idx.cut,]
    temp.ind = temp.ind+1
  }
  return(filter_result)
}





#all snp.result
all.result <- rbind(known.result[idx.fil,],dis.result)
#load in the heritability estimate
heri.est <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/BCAC_heritability.csv")
#reorder heri.est 
heri.est <- heri.est[c(2,4,5,3,1),c(2,4,5,3,1)]
gwas.heri <- diag(as.matrix(heri.est))

#five subtypes
all.snp.explain <- rep(0,5)

for(i in 1:5){
  all.snp.explain[i] = SigmaEst(all.result[,2*i-1],all.result[,2*i],all.result[,13]) 
}
############temp results for the Lumianl A like and standard log odds ratio of all the 210 identified SNPs
write.csv(all.result,file="./heritability/result/all_snps_log_odds.csv")


#overall
i <- 6
all.snp.explain.overall <- SigmaEst(all.result[,2*i-1],all.result[,2*i],all.result[,13]) 
print(all.snp.explain.overall)
new.snp.explain <- rep(0,5)

for(i in 1:5){
  new.snp.explain[i] = SigmaEst(dis.result[,2*i-1],dis.result[,2*i],dis.result[,13]) 
}
#overall
i <- 6
new.snp.explain.overall <- SigmaEst(dis.result[,2*i-1],dis.result[,2*i],dis.result[,13]) 
print(new.snp.explain.overall)
#the heritability explained for two-fold increase of family risk

#overall heritability based on oncoarray paper
#summarize the result
result.sum <- cbind(c(all.snp.explain.overall,all.snp.explain),
                    c(new.snp.explain.overall,new.snp.explain),
                    c(2*log(2)*0.41,gwas.heri),
                    c(all.snp.explain.overall,all.snp.explain)/c(2*log(2)*0.41,gwas.heri))
row.names(result.sum) <- c("overall",
                           "Luminial A","Luminal B",
                           "Luminal B HER2-",
                           "HER2 Enriched",
                           "Triple Negative")
colnames(result.sum) <- c("simga_all_identified_snps",
                          "sigma_new_identified_snps",
                          "sigma_gwas",
                          "all_gwas_proportion")
write.csv(result.sum,file = "./heritability/result/result.sum.csv")

############PRS for the top 1% increase risk compared to the standard
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
#icog.train.id <- split.id[[1]]
#onco.train.id <- split.id[[2]]
onco.test.id <- split.id[[5]]



data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1,data2$Grade1)
onco.test <- c(1:nrow(y.pheno.mis2))
#onco.test <- which(data2[,1]%in%onco.test.id)
y.pheno.mis2 <- y.pheno.mis2[onco.test,]
#y.pheno.mis2 <- cbind(data2$Behaviour1,data2$PR_status1,data2$ER_status1,data2$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour","PR",
                           "ER","HER2","Grade")

x.test.all.mis2_known <- data2[onco.test,c(27:203,205)]
########all the 210 SNPs
x.test.all.mis_dis <- x.test.all.mis.dis[onco.test,]
x.test.all.mis2.all <- as.matrix(cbind(x.test.all.mis2_known,x.test.all.mis_dis))
n <- nrow(y.pheno.mis2)
#####PRS for everyone in OncoArray
PRS <- matrix(0,n,6)
#PRS_temp <- PRS[,1]
for(i in 1:6){
  PRS[,i] <- x.test.all.mis2.all%*%all.result[,2*i-1]
}
#####risk of top 1% vs. population average 

RelativeRisk <- function(PRS_temp){
  idx.control <- which(y.pheno.mis2[,1]==0)
  PRS_temp_control <- PRS_temp[idx.control]
  idx <- which(PRS_temp >= quantile(PRS_temp_control,c(0.48))& PRS_temp <= quantile(PRS_temp_control,c(0.52)))
  odds1 <- sum(y.pheno.mis2[idx,1])/(length(idx)-sum(y.pheno.mis2[idx,1]))
  idx2 <- which(PRS_temp >= quantile(PRS_temp,0.99))
  odds2 <- sum(y.pheno.mis2[idx2,1])/(length(idx2)-sum(y.pheno.mis2[idx2,1]))
  or <- odds2/odds1
  p <- 0.124
  rr <- or/(1-p+(p*or))
  return(rr)
}
idx.case <- which(y.pheno.mis2[,1]==1)
PRS_temp <- PRS[,6]
idx.control <- which(y.pheno.mis2[,1]==0)
PRS_temp_case <- PRS_temp[idx.case]
PRS_temp_control <- PRS_temp[idx.control]
qnorm(0.99,mean(PRS_temp_control),sd(PRS_temp_control))
qnorm(0.01,mean(PRS_temp_control),sd(PRS_temp_control))
result <- pnorm(1.376025,mean=mean(PRS_temp_control)+var(PRS_temp_control),sd = sd(PRS_temp_control),lower.tail = F)/pnorm(1.376025,mean=mean(PRS_temp_control),sd = sd(PRS_temp_control),lower.tail = F)

result2 <- pnorm(-1.437285,mean=mean(PRS_temp_control)+var(PRS_temp_control),sd = sd(PRS_temp_control),lower.tail = T)/pnorm(-1.437285,mean=mean(PRS_temp_control),sd = sd(PRS_temp_control),lower.tail = T)



# 
# mean(PRS_temp_case)
# sd(PRS_temp_case)
relative.risk.all <- rep(0,6)
for(i in 1:6){
  relative.risk.all[i] <- RelativeRisk(PRS[,i])
}
