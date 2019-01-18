#-------------------------------------------------------------------
# Update Date: 01/18/2019
# Create Date: 01/18/2019
# Goal: estimate the log odds ratio and variance for 32 discovery SNPs
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
result <- NULL
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')

for(i1 in 1:35){
  load(paste0("./discovery_SNP/additive_model/result/intrinsic_subtype_",i1,".Rdata"))
  result <- rbind(result,test.result.second.wald)
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

###remove 3 SNP after conditional p-value threshold becomes 1e-06
result <- result[-c(7,24,33),]
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




















