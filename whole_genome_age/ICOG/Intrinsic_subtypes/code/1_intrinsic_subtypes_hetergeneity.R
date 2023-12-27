args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
#load the data
load(paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_082119.Rdata"))
startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}
size = 1000
num = nrow(meta_result_shared_1p)
start_end_result = startend(num,size,i1)
start = start_end_result[1]
end = start_end_result[2]
C = matrix(c(rep(1,4),
             c(-1,0,0,0),
             c(0,-1,0,0),
             c(0,0,-1,0),
             c(0,0,0,-1)),4,5)
p_sub = rep(0,end-start+1)
temp = 1
for(k in start:end){
  beta = as.numeric(meta_result_shared_1p[k,16:20])
  beta_cov = matrix(as.numeric(meta_result_shared_1p[k,21:45]),5,5)
  beta_heter = C%*%beta
  cov_heter = C%*%beta_cov%*%t(C)
  chi_test = t(beta_heter)%*%solve(cov_heter)%*%beta_heter
  p_sub[temp] = pchisq(chi_test,length(beta_heter),lower.tail = T)
  temp = temp + 1
}
save(p_sub, file = paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/p_heter_sub_",i1,".rdata"))
