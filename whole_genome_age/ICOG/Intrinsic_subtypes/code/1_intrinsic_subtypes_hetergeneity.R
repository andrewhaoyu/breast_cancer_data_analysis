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
C2 = matrix(c(1,-1),1,2)
p_sub = rep(0,end-start+1)
p_sub2 = rep(0,end-start+1)
temp = 1

idx <- which(meta_result_shared_1p$var_name=="19_17354586_C_A")

for(k in start:end){
  beta = as.numeric(meta_result_shared_1p[k,16:20])
  beta_cov = matrix(as.numeric(meta_result_shared_1p[k,21:45]),5,5)
  beta_heter = C%*%beta
  cov_heter = C%*%beta_cov%*%t(C)
  chi_test = t(beta_heter)%*%solve(cov_heter)%*%beta_heter
  p_sub[temp] = pchisq(chi_test,length(beta_heter), lower.tail = F)
  
  beta_sub = beta[c(1,5)]
  beta_cov_sub = beta_cov[c(1,5),c(1,5)]
  beta_heter = C2%*%beta_sub
  cov_heter = C2%*%beta_cov_sub%*%t(C2)
  chi_test = t(beta_heter)%*%solve(cov_heter)%*%beta_heter
  p_sub2[temp] = pchisq(chi_test,length(beta_heter), lower.tail = F)
  temp = temp + 1
}
p_list = list(p_sub, p_sub2)
save(p_list, file = paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/p_heter_sub_",i1,".rdata"))
