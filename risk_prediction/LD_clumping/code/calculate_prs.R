#-------------------------------------------------------------------
# Update Date: 11/22/2018
# Create Date: 11/22/2018
# Goal: use plink to calculate prs
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
prs.code <- "/spin1/users/zhangh24/plink --score /spin1/users/zhangh24/BCAC/prs_file/Lu_standard_prs.file no-sum no-mean-imputation --dosage /spin1/users/zhangh24/BCAC/impute_onco_dosage/dosage_all noheader skip0=1 skip1=1 format=1 --fam /spin1/users/zhangh24/BCAC/impute_onco/onco_plink.fam --out /spin1/users/zhangh24/BCAC/prs_out/Lu_standard_prs7"


write.table(prs.code,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/caculate.prs.sh"),col.names = F,row.names = F,quote=F)
