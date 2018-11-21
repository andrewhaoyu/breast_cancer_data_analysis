#-------------------------------------------------------------------
# Update Date: 11/21/2018
# Create Date: 11/21/2018
# Goal: convert chromsome filt to PRS
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------

convert_dosage <- rep("c",22)
convert_dosage <- data.frame(convert_dosage,stringsAsFactors=F)
for(i in 1:22){
  convert_dosage_code <- paste0("zcat /spin1/users/zhangh24/BCAC/impute_onco/chr$i.gz | awk -v C=2 -v TF=/spin1/users/zhangh24/BCAC/impute_plink_onco/clump_snp -f /home/zhangh24/bin/Element.Matching.awk  | awk -v chr=",i," '{printf chr\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5; for(i=6; i<NF; i+=3) {if($(i+0) == 0 && $(i+1) == 0 && $(i+2) == 0) printf \"\\tNA\"; else printf \"\\t\"$(i+0)*2+$(i+1)*1+$(i+2)*0}; printf \"\\n\"}' | gzip > /spin1/users/zhangh24/BCAC/impute_onco_dosage/dosage_chr",i)
  convert_dosage[i,1] <- convert_dosage_code
}
#write out the command and submit use cluster
write.table(convert_dosage,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/convert.dosage.sh"),col.names = F,row.names = F,quote=F)

