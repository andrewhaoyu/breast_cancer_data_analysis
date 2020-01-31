#-------------------------------------------------------------------
# Update Date: 01/22/2020
# Create Date: 11/20/2018
# Goal: LD pruning
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
# r2thr = 0.05; kbpthr = 1000; pthr=0.1
# #trait_name = "BCAC_metacase"
# #---------------------------------------#---------------------------------------
# LD.clump <- rep("c",22)
# LD.clump <- data.frame(LD.clump,stringsAsFactors=F)
# for(i in 1:22){
#   LD.clump.code <- paste0("/data/zhangh24/plink --bfile /data/zhangh24/BCAC/impute_plink_onco/chr",i,"_plink --clump /data/zhangh24/BCAC/impute_plink_onco/LD_assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump")
#   LD.clump[i,1] <- LD.clump.code
# }
# #write out the command and submit use cluster
# write.table(LD.clump,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/LD.clump.sh"),col.names = F,row.names = F,quote=F)
# 
# 
# 
# 
# #merge the LD clumped SNP list
# 
# total <- 0
# for(i in 1:22){
#   data <- read.table(paste0("/data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump.clumped"),header=T)
#   temp <- nrow(data)
#   total <- temp+total
# }
# 
# SNP <- rep("c",total)
# total <- 0
# for(i in 1:22){
#   data <- read.table(paste0("/data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump.clumped"),header=T)
#   temp <- nrow(data)
#   SNP[total+(1:temp)] <- as.character(data[,3])
#   total <- total+temp
# }
# write.table(SNP,file = "/data/zhangh24/BCAC/impute_plink_onco/clump_snp",
#             row.names = F,col.names = F,quote=F)


#data <- fread("/data/zhangh24/BCAC/impute_plink_onco/clump_snp")


# system('rm /data/zhangh24/BCAC/impute_plink_onco/clump_snp')
# #create a file called clump_snp using touch command
# system('touch /data/zhangh24/BCAC/impute_plink_onco/clump_snp')
# #take out the SNPs using awk command
# system("awk '{for(i=1;i<=22;i++)
# {print $3} /data/zhangh24/BCAC/impute_plink_onco/chr$i_ld_clump.clumped >> /data/zhangh24/BCAC/impute_plink_onco/clump_snp}'
#        done")
 




r2thr = 0.1; kbpthr = 500; pthr=0.1
#trait_name = "BCAC_metacase"
#---------------------------------------#---------------------------------------
LD.clump <- rep("c",22)
LD.clump <- data.frame(LD.clump,stringsAsFactors=F)
for(i in 1:22){
  LD.clump.code <- paste0("/data/zhangh24/plink --bfile /data/zhangh24/BCAC/impute_plink_onco/chr",i,"_plink --clump /data/zhangh24/BCAC/impute_plink_onco/LD_assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump_121019")
  LD.clump[i,1] <- LD.clump.code
}
#write out the command and submit use cluster
write.table(LD.clump,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/LD.clump_121019.sh"),col.names = F,row.names = F,quote=F)




#merge the LD clumped SNP list

total <- 0
for(i in 1:22){
  data <- read.table(paste0("/data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump_121019.clumped"),header=T)
  temp <- nrow(data)
  total <- temp+total
}

SNP <- rep("c",total)
CHR <- rep(0,total)
BP <- rep(0,total)
total <- 0
for(i in 1:22){
  data <- read.table(paste0("/data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump_121019.clumped"),header=T)
  temp <- nrow(data)
  SNP[total+(1:temp)] <- as.character(data[,3])
  CHR[total+(1:temp)] <- as.numeric(data[,1])
  BP <- as.numeric(data[,4])
  total <- total+temp
}

#check whether the 313 SNPs are within the LD clumping set
setwd('/data/zhangh24/breast_cancer_data_analysis/')
load("./data/Nasim_313SNPs_complete_information.Rdata")
head(snp.new)
#find the SNPs that are dropped by the LD-pruning procedure
idx <- which(snp.new$SNP.ONCO%in%
               SNP!=T)
length(idx)
library(data.table)
assoc <- as.data.frame(fread("/data/zhangh24/BCAC/impute_plink_onco/LD_assoc"))

for(k in 1:length(idx)){
  #find all the 313 SNPs that were filttered due to LD-pruning procedure
  temp.SNP <- snp.new$SNP.ONCO[idx[k]]
  #find all the SNPs +-500kb in it and remove them
  zdx <- which(assoc$SNP==temp.SNP)
  
  ldx <- which(assoc$CHR==assoc$CHR[zdx]&
                 assoc$BP>=assoc$BP[zdx]-500000&
                 assoc$BP<=assoc$BP[zdx]+500000&assoc$SNP%in%snp.new$SNP.ONCO==F)
  if(length(ldx)!=0){
    print(k)
    remove.list <- assoc$SNP[ldx]
    cdx <- which(SNP%in%remove.list)
    SNP <- SNP[-cdx]
  }
  SNP <- c(SNP,temp.SNP)
}
#check whether all the 313 SNPs are in the list
idx <- which(snp.new$SNP.ONCO%in%
               SNP!=T)
length(idx)

write.table(SNP,file = "/data/zhangh24/BCAC/impute_plink_onco/clump_snp_121019",
            row.names = F,col.names = F,quote=F)

#data <- fread("/data/zhangh24/BCAC/impute_plink_onco/clump_snp")


# system('rm /data/zhangh24/BCAC/impute_plink_onco/clump_snp')
# #create a file called clump_snp using touch command
# system('touch /data/zhangh24/BCAC/impute_plink_onco/clump_snp')
# #take out the SNPs using awk command
# system("awk '{for(i=1;i<=22;i++)
# {print $3} /data/zhangh24/BCAC/impute_plink_onco/chr$i_ld_clump.clumped >> /data/zhangh24/BCAC/impute_plink_onco/clump_snp}'
#        done")















