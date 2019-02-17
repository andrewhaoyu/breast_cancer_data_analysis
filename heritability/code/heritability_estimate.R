#-------------------------------------------------------------------
# Update Date: 01/20/2019
# Create Date: 01/18/2019
# Goal: estimate heritability for breast cancer overall risk and subtypes risk
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
Outvec <- function(heter.result){
  logodds <- heter.result[[1]]
  var.log <- diag(heter.result[[2]])
  standard.log <- heter.result[[3]]
  var.standard <- heter.result[[4]]
  result.vec <- rep(0,13)
  temp = 1
  for(i in 1:5){
    result.vec[temp] = logodds[i]
    temp = temp+1
    result.vec[temp] = var.log[i]
    temp = temp+1
  }
  result.vec[temp] = standard.log
  temp = temp + 1
  result.vec[temp] = var.standard
  temp = temp + 1
  result.vec[temp] = heter.result[[5]]
  return(result.vec)
}


result <- NULL
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
result_standard <- matrix(0,35,2)
result <- matrix(0,35,10)
freq = rep(0,35)
for(i1 in 1:35){
  load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_heter_herita",i1,".Rdata"))
 # load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_herita_",i1,".Rdata"))
  temp = Outvec(test.result.second.wald)
  result[i1,] <- temp[1:10]
  freq[i1] <- temp[13]
  #######plug in the log odds ratio and var from standard logistic regression
  result_standard[i1,1] <- temp[11]
  result_standard[i1,2] <- temp[12]
  
}
#SNP <- c(colnames(icog.julie),colnames(discovery.snp.icog)[1:18])

##################discovery snp were ordered based on the order they are extracted
discovery_snp <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv",header=T)
#SNP <- discovery_snp$SNP.ICOGS

################match the discovery snps to the order in the paper
discovery_snp_paper_order <- read.csv("./data/discovery_snp_paper_order.csv",header=T)
chr.pos.paper <- paste0(discovery_snp_paper_order$CHR,":",discovery_snp_paper_order$position)


chr.pos <- paste0(discovery_snp$CHR.x,":",discovery_snp$position)

idx.match <- match(chr.pos.paper,
                   chr.pos)
discovery_snp_new <- discovery_snp[idx.match,]
result <- result[idx.match,]
result_standard <- result_standard[idx.match,]
freq = freq[idx.match]
###remove 3 SNP after conditional p-value threshold becomes 1e-06
discovery_snp_new <- discovery_snp_new[-c(7,24,33),]
result <- result[-c(7,24,33),]
result_standard <- result_standard[-c(7,24,33),]
freq = freq[-c(7,24,33)]
#load the results from Oncoarray paper
all_result <- fread("./data/oncoarray_bcac_public_release_oct17.txt",header=T)
all_result <- all_result %>% 
  mutate(chrpos=paste0(chr,":",position_b37))
discovery_snp_new <- left_join(discovery_snp_new,all_result,
                               by= "var_name")
# result_standard <- cbind(as.numeric(discovery_snp_new$bcac_onco_icogs_gwas_beta),
#                          as.numeric(discovery_snp_new$bcac_onco_icogs_gwas_se^2))
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
  mutate(chrpos = paste0(CHR,":",position))
fine_map_new <- left_join(fine_map,all_result,by="chrpos")
#one snp duplicated
fine_map_new <- fine_map_new[-c(62,120),]
known.result[,c(11,12)] <- cbind(as.numeric(fine_map_new$bcac_onco2_beta),as.numeric(fine_map_new$bcac_onco2_se)^2)
known.result[178,c(11,12)] <- temp.result



#all snp.result
all.result <- rbind(known.result,dis.result)
#load in the heritability estimate
heri.est <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/BCAC_heritability.csv")
#reorder heri.est 
heri.est <- heri.est[c(2,4,5,3,1),c(2,4,5,3,1)]
gwas.heri <- diag(as.matrix(heri.est))
SigmaEst <- function(beta,sigma,p){
  n <- length(p)
  result <- sum(2*p*(1-p)*(beta^2-sigma))
  return(result)
}
#five subtypes
all.snp.explain <- rep(0,5)

for(i in 1:5){
  all.snp.explain[i] = SigmaEst(all.result[,2*i-1],all.result[,2*i],all.result[,13]) 
}

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
