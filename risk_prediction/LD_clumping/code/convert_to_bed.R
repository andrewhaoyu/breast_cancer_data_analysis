#-------------------------------------------------------------------
# Update Date: 11/14/2018
# Create Date: 11/14/2018
# Goal: convert impute gen fortmat to bed format
# Author: Haoyu Zhang
#-------------------------------------------------------------------

convert.chr <- rep("c",22)
convert.chr <- data.frame(convert.chr,stringsAsFactors=F)

for(i in 1:22){
  convert.code <- paste0("/spin1/users/zhangh24/plink --gen /spin1/users/zhangh24/BCAC/impute_onco/chr",i,".gz --sample /spin1/users/zhangh24/BCAC/sample_onco.txt --oxford-single-chr ",i," -make-bed --out /spin1/users/zhangh24/BCAC/impute_plink_onco/chr",i,"_plink")
  convert.chr[i,1] <- convert.code
}


#write out the command and submit use cluster
write.table(convert.chr,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/convert.chr.sh"),col.names = F,row.names = F,quote=F)




