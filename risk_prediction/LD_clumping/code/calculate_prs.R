#-------------------------------------------------------------------
# Update Date: 11/24/2018
# Create Date: 11/22/2018
# Goal: use plink to calculate prs for different subtyeps
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
subtypes <- c("Luminial_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
select.names <- c(subtypes,paste0("eb_",subtypes))
select.names <- c("standard",
                  select.names)
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
            1E-04,5E-04,1E-03,5E-03,1E-02)
#create the file for standard logisic regression
n.pthres <- length(pthres)


prs.code <- rep("c",length(select.names))
prs.code <- data.frame(prs.code,stringsAsFactors=F)
temp <- 1
for(i in 1:n.pthres){
  for(j in 1:length(select.names)){
    prs.code.temp <- paste0("/spin1/users/zhangh24/plink --score /spin1/users/zhangh24/BCAC/prs_file/",select.names[j],"_prs_",i,".file no-sum no-mean-imputation --dosage /spin1/users/zhangh24/BCAC/impute_onco_dosage/dosage_all noheader skip0=1 skip1=1 format=1 --fam /spin1/users/zhangh24/BCAC/impute_onco/onco_plink.fam --out /spin1/users/zhangh24/BCAC/prs_out/",select.names[j],"_prs_",i,"_out")
    prs.code[temp,1] <- prs.code.temp
    temp <- temp+1
  }
  }
#write out the command and submit use cluster
write.table(prs.code,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/caculate.prs.sh"),col.names = F,row.names = F,quote=F)








