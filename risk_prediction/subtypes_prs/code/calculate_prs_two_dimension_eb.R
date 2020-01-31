#-------------------------------------------------------------------
# Update Date: 11/24/2018
# Create Date: 11/22/2018
# Goal: use plink to calculate prs for different subtyeps
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
# subtypes <- c("Luminal_A",
#               "Luminal_B",
#               "Luminal_B_HER2Neg",
#               "HER2_Enriched",
#               "TN")
# #select.names <- c(subtypes,paste0("eb_",subtypes))
# #select.names <- c("standard",subtypes)
# select.names <- c(subtypes)
# # pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
# #             1E-04,5E-04,1E-03,5E-03,1E-02)
# pthres <- c(1E-30,1E-10,5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02)
# 
# #create the file for standard logisic regression
# n.pthres <- length(pthres)
# 
# 
# prs.code <- rep("c",length(select.names))
# prs.code <- data.frame(prs.code,stringsAsFactors=F)
# temp <- 1
# for(i in 1:n.pthres){
#   for(k in 1:n.pthres){
#   for(j in 1:length(select.names)){
#     # noheader means the dosage file has no header
#     # skip0=N      Number of fields to skip before SNP
#     # skip1=N       Number of fields to skip between SNP and A1
#     # format=N      Dosage, two probabilities or three (N=1,2,3)
#     prs.code.temp <- paste0("/data/zhangh24/plink --score /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,"_",k,".file no-sum no-mean-imputation --dosage /data/zhangh24/BCAC/impute_plink_onco/dosage_all_new_try noheader skip0=1 skip1=1 format=1 --fam /data/zhangh24/BCAC/impute_onco/onco_plink.fam --out /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_",i,"_",k,"_out")
#     prs.code[temp,1] <- prs.code.temp
#     temp <- temp+1
#   }
#   }
# }
# #write out the command and submit use cluster
# write.table(prs.code,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/code/caculate.prs.2d.sh"),col.names = F,row.names = F,quote=F)







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
# pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
#             1E-04,5E-04,1E-03,5E-03,1E-02)
pthres <- c(1E-30,1E-10,5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02)

#create the file for standard logisic regression
n.pthres <- length(pthres)


prs.code <- rep("c",length(select.names))
prs.code <- data.frame(prs.code,stringsAsFactors=F)
temp <- 1
for(i in 1:n.pthres){
  for(k in 1:n.pthres){
  for(j in 1:length(select.names)){
    # noheader means the dosage file has no header
    # skip0=N      Number of fields to skip before SNP
    # skip1=N       Number of fields to skip between SNP and A1
    # format=N      Dosage, two probabilities or three (N=1,2,3)
    prs.code.temp <- paste0("/data/zhangh24/plink --score /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_pvaluecut_",i,"_",k,"_121019_eb.file no-sum no-mean-imputation --dosage /data/zhangh24/BCAC/impute_plink_onco/dosage_all_new_try_121019 noheader skip0=1 skip1=1 format=1 --fam /data/zhangh24/BCAC/impute_onco/onco_plink.fam --out /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/",select.names[j],"_prs_",i,"_",k,"_out_121019_eb")
    prs.code[temp,1] <- prs.code.temp
    temp <- temp+1
  }
  }
}
#write out the command and submit use cluster
write.table(prs.code,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/code/caculate.prs.2d.sh_121019_eb"),col.names = F,row.names = F,quote=F)







