#-------------------------------------------------------------------
# Update Date: 11/17/2018
# Create Date: 11/16/2018
# Goal: convert sample data into BGEN
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#prepare sample for all of the subjects
names.subtypes <-  c("Luminal_A","Luminal_B",
                     "Luminal_B_HER2Neg",
                     "HER2Enriched",
                     "TripleNeg")

convert.to.bgen <- rep("c",length(names.subtypes))
convert.to.bgen <- data.frame(convert.to.bgen,stringsAsFactors=F)

for(i in 1:length(names.subtypes)){
  convert.to.bgen.code <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g /spin1/users/zhangh24/BCAC/impute_onco/onco_all.gz -s /spin1/users/zhangh24/BCAC/sample_onco.txt -incl-samples /spin1/users/zhangh24/BCAC/test_sample_",names.subtypes[i],".txt -incl-rsids /spin1/users/zhangh24/BCAC/impute_onco/onco_1p_shared_id.txt -og /spin1/users/zhangh24/BCAC/impute_test_bgen/",names.subtypes[i],"_test.bgen")
  convert.to.bgen[i,1] <- convert.to.bgen.code
}


write.table(convert.to.bgen,file="/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/convert_to_bgen",
            row.names=F,col.names = F,
            quote=F)







