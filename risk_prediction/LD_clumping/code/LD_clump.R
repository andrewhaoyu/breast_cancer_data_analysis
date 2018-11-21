#-------------------------------------------------------------------
# Update Date: 11/21/2018
# Create Date: 11/20/2018
# Goal: LD pruning
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------
r2thr = 0.1; kbpthr = 250; pthr=0.1
#trait_name = "BCAC_metacase"
#---------------------------------------#---------------------------------------
LD.clump <- rep("c",22)
LD.clump <- data.frame(LD.clump,stringsAsFactors=F)
for(i in 1:22){
  LD.clump.code <- paste0("/spin1/users/zhangh24/plink --bfile /spin1/users/zhangh24/BCAC/impute_plink_onco/chr",i,"_plink --clump /spin1/users/zhangh24/BCAC/impute_plink_onco/LD_assoc --clump-p1 ",pthr," --clump-p2 ", pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /spin1/users/zhangh24/BCAC/impute_plink_onco/chr",i,"_ld_clump")
  LD.clump[i,1] <- LD.clump.code
}
#write out the command and submit use cluster
write.table(LD.clump,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/LD.clump.sh"),col.names = F,row.names = F,quote=F)




#merge the LD clumped SNP list
system('rm /spin1/users/zhangh24/BCAC/impute_plink_onco/clump_snp')
system('touch /spin1/users/zhangh24/BCAC/impute_plink_onco/clump_snp')
system("awk '{print $3}' /spin1/users/zhangh24/BCAC/impute_plink_onco/chr$i_ld_clump.clumped >> /spin1/users/zhangh24/BCAC/impute_plink_onco/clump_snp
       done")
 









