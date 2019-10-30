# args = commandArgs(trailingOnly = T)
# i1 = as.numeric(args[[1]])


load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_shared_1p_sex.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_shared_1p_sex.Rdata")
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
second.num <- 5


icog_score_infor <- icog_result_shared_1p[,11:(11+second.num+second.num^2-1)]
onco_score_infor <- onco_result_shared_1p[,11:(11+second.num+second.num^2-1)]



icog_onco_score_infor <- cbind(icog_score_infor,onco_score_infor)



n <- nrow(icog_onco_score_infor)




#registerDoParallel(no.cores)
#pvalue <- rep(0,nrow(icog_onco_score_infor))


#pvalue.list <- foreach(i=1:size)%dopar%
#{
#print(i)
start <- 1
end <- nrow(icog_onco_score_infor)
# start.end <- startend(n,size,i1)
# start <- start.end[1]
# end <- start.end[2]
pvalue_sub <- rep(0,end-start+1)
temp = 1
for(j in start:end){
  print(j)
  
  icog_onco_score_infor_oneline <- icog_onco_score_infor[j,]
  pvalue_sub[temp] <- MetaPfunction(icog_onco_score_infor_oneline,second.num)
  temp = temp+1
}
#since sex chromosome is small
#just directly run it interactively
p.value <- pvalue_sub

load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_result_shared_1p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_result_shared_1p.Rdata")


meta_result_shared_1p_no_pvalue <- icog_result_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]



meta_result_shared_1p <- cbind(meta_result_shared_1p_no_pvalue,p.value)
save(meta_result_shared_1p,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_sex.Rdata"))
meta_result_shared_1p_sex <- meta_result_shared_1p

load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata"))

head(meta_result_shared_1p)
#match the column name between the sex dataset and #22 auto chromosomes data
names.sex <- colnames(meta_result_shared_1p_sex)
names.standard <- colnames(meta_result_shared_1p)
#sex data 1:3,6:15
#standard data 1:8,11:15
meta_result_shared_1p_final = rbind(
  meta_result_shared_1p[,c(1:8,11:15)],
  meta_result_shared_1p_sex[,c(1:3,6:15)] )
meta_result_shared_1p <- meta_result_shared_1p_final
save(meta_result_shared_1p,file = "/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_final.Rdata")
meta_result_shared_1p <- meta_result_shared_1p[,c(2,9,3,13)]
write.table(meta_result_shared_1p,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_fixed.txt"),row.names = F,
            quote=F)


idx <- which(meta_result_shared_1p$CHR==12&
               meta_result_shared_1p$position==111600134)
