#-------------------------------------------------------------------
# Update Date: 12/04/2018
# Create Date: 11/20/2018
# Goal: LD pruning
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
r2thr = 0.05; kbpthr = 1000; pthr=0.1
#trait_name = "BCAC_metacase"
#---------------------------------------#---------------------------------------
LD.clump <- rep("c",22)
LD.clump <- data.frame(LD.clump,stringsAsFactors=F)
for(i in 1:22){
  LD.clump.code <- paste0("/data/zhangh24/plink --bfile /data/zhangh24/BCAC/impute_plink_onco/chr",i,"_plink --clump /data/zhangh24/BCAC/impute_plink_onco/LD_assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump")
  LD.clump[i,1] <- LD.clump.code
}
#write out the command and submit use cluster
write.table(LD.clump,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/LD.clump.sh"),col.names = F,row.names = F,quote=F)




#merge the LD clumped SNP list

total <- 0
for(i in 1:22){
  data <- read.table(paste0("/data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump.clumped"),header=T)
  temp <- nrow(data)
  total <- temp+total
}

SNP <- rep("c",total)
total <- 0
for(i in 1:22){
  data <- read.table(paste0("/data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump.clumped"),header=T)
  temp <- nrow(data)
  SNP[total+(1:temp)] <- as.character(data[,3])
  total <- total+temp
}
write.table(SNP,file = "/data/zhangh24/BCAC/impute_plink_onco/clump_snp",
            row.names = F,col.names = F,quote=F)


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
total <- 0
for(i in 1:22){
  data <- read.table(paste0("/data/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump_121019.clumped"),header=T)
  temp <- nrow(data)
  SNP[total+(1:temp)] <- as.character(data[,3])
  total <- total+temp
}
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















