#-------------------------------------------------------------------
# Update Date: 11/24/2018
# Create Date: 11/22/2018
# Goal: use plink to calculate prs for different subtyeps
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
subtypes <- c("Luminal_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
#select.names <- c(subtypes,paste0("eb_",subtypes))
#select.names <- c("standard",subtypes)
select.names <- c(subtypes)
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
            1E-04,5E-04,1E-03,5E-03,1E-02)
#create the file for standard logisic regression
n.pthres <- length(pthres)


prs.code <- rep("c",length(select.names))
prs.code <- data.frame(prs.code,stringsAsFactors=F)
temp <- 1
for(i in 1:n.pthres){
  for(j in 1:length(select.names)){
    # noheader means the dosage file has no header
    # skip0=N      Number of fields to skip before SNP
    # skip1=N       Number of fields to skip between SNP and A1
    # format=N      Dosage, two probabilities or three (N=1,2,3)
    prs.code.temp <- paste0("/data/zhangh24/plink --score /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,".file no-sum no-mean-imputation --dosage /data/zhangh24/BCAC/impute_onco_dosage/dosage_all noheader skip0=1 skip1=1 format=1 --fam /data/zhangh24/BCAC/impute_onco/onco_plink.fam --out /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_",i,"_out")
    prs.code[temp,1] <- prs.code.temp
    temp <- temp+1
  }
  }
#write out the command and submit use cluster
write.table(prs.code,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/code/caculate.prs.sh"),col.names = F,row.names = F,quote=F)


# temp = read.table("/data/zhangh24/BCAC/impute_onco/onco_plink.fam",header=F)
# temp2 = read.table("/data/zhangh24/test/sample.txt",
#                    header=T)
# idx <- which(as.character(temp$V1)!=as.character(temp2$ID_1)[2:nrow(temp2)])


#"/data/zhangh24/plink --score /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/test.file no-sum no-mean-imputation --map /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/onco_map.txt --dosage /data/zhangh24/BCAC/impute_onco_dosage/dosage_all noheader skip0=1 skip1=1 format=1 --fam /data/zhangh24/BCAC/impute_onco/onco_plink.fam --out /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/test_out"

#"/data/zhangh24/plink --map /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/onco_map.txt --dosage /data/zhangh24/BCAC/impute_onco_dosage/dosage_all noheader skip0=1 skip1=1 format=1 --fam /data/zhangh24/BCAC/impute_onco/onco_plink.fam --out /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/dosage_test_out"

idx <- which(onco_info$rs_id=="rs370540207:121450795:T:G")