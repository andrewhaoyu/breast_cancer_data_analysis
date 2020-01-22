#-------------------------------------------------------------------
# Update Date: 11/14/2018
# Create Date: 11/14/2018
# Goal: merge the imputed gen format data of BCAC as one single data
# Author: Haoyu Zhang
#-------------------------------------------------------------------
merge.code <- paste0("cat ")
for(i in 1:22){
  #take out all the files with chri with grep command
  #fixed option = T only take out the exact match
  
  
    merge.code <- paste0(" ",
                         merge.code,
                         "/data/zhangh24/BCAC/impute_onco/chr",
                         i,
                         ".gz",
                         " ")  
  


  
}
merge.code <- paste0(merge.code,
                     "> /data/zhangh24/BCAC/impute_onco/onco_all",
                     ".gz")
system(merge.code)
 #write out the command and submit use cluster
write.table(merge.code,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/merge.gen.to.one.onco.sh"),col.names = F,row.names = F,quote=F)
# 



merge.code <- paste0("cat ")
for(i in 1:22){
  #take out all the files with chri with grep command
  #fixed option = T only take out the exact match
  
  
  merge.code <- paste0(" ",
                       merge.code,
                       "/data/zhangh24/BCAC/impute_icogs/chr",
                       i,
                       ".gz",
                       " ")  
  
  
  
  
}
merge.code <- paste0(merge.code,
                     "> /data/zhangh24/BCAC/impute_icogs/icogs_all_",
                     "gen.gz")
#system(merge.code)
#write out the command and submit use cluster
write.table(merge.code,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/LD_clumping/code/merge.gen.to.one.icogs.sh"),col.names = F,row.names = F,quote=F)

