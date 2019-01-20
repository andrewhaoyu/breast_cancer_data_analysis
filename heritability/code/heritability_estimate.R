#-------------------------------------------------------------------
# Update Date: 01/20/2019
# Create Date: 01/18/2019
# Goal: estimate heritability for breast cancer overall risk and subtypes risk
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
result <- NULL
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
result_standard <- matrix(0,35,2)
for(i1 in 1:35){
  load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_heter_herita",i1,".Rdata"))
 # load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_herita_",i1,".Rdata"))
  result <- rbind(result,test.result.second.wald[[1]])
  #######plug in the log odds ratio and var from standard logistic regression
  result_standard[i1,1] <- test.result.second.wald[[2]]
  result_standard[i1,2] <- test.result.second.wald[[3]]
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
###remove 3 SNP after conditional p-value threshold becomes 1e-06
result <- result[-c(7,24,33),]
result_standard <- result_standard[-c(7,24,33),]
###convert the result to log odds ratio
n.snp <- 32
n.subtypes <- 5
dis.result <- matrix(0,n.snp,2*n.subtypes)
for(i in 1:n.snp){
  for(j in c(2*(1:n.subtypes)-1)){
    dis.result[i,j] <- log(as.numeric(strsplit(result[i,j],"\\(")[[1]][1]))
    
  }
}
dis.result[,c(2*(1:n.subtypes))] <-as.matrix(result[,c(2*(1:n.subtypes))])
dis.result <- as.data.frame(dis.result)

colnames(dis.result)[c(2*(1:n.subtypes))-1] <- paste0("logodds_",c("Luminial A","Luminal B",
                        "Luminal B HER2-",
                        "HER2 Enriched",
                        "Triple Negative"))
colnames(dis.result)[c(2*(1:n.subtypes))] <- paste0("var_",c("Luminial A","Luminal B",
                                                                   "Luminal B HER2-",
                                                                   "HER2 Enriched",
                                                                   "Triple Negative"))
colnames(result_standard) <- c("logoddsd_standard",
                               "var_standard")
dis.result <- cbind(dis.result,result_standard)

##load the results for 178 known SNPs
known.result <- matrix(0,178,12)
colnames(known.result) <- colnames(dis.result)

Outvec <- function(heter.result){
  logodds <- heter.result[[1]]
  var.log <- diag(heter.result[[2]])
  standard.log <- heter.result[[3]]
  var.standard <- heter.result[[4]]
  result.vec <- rep(0,12)
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
  return(result.vec)
}


for(i in 1:178){
  print(i)
  load(paste0("./known_SNPs/known_SNPs_analysis_G_revised/intrinsic_subtypes_pc_additive/result/heter_result_origin",i,".Rdata"))
  known.result[i,] <- Outvec(heter.result)
}

#load in the heritability estimate
















