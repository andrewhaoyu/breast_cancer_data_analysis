#meta analysis for two interaction model
arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
load("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_fixed_5p_two_interaction.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_random_5p_two_interaction.Rdata")
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")

num = nrow(onco_result_fixed_5p)
size = 1000
start.end <- startend(num,size,i1)
start <- start.end[1]
end <- start.end[2]
file.num <- end-start+1
p_sub1 <- rep(0,file.num)
p_sub2 <- rep(0,file.num)
temp = 1
n.support = 2
n.test = 9
n.fixed = n.support+n.test
for(i in start:end){
  score = t(as.numeric(onco_result_fixed_5p[i,10+(1:n.fixed)]))
  
  infor = matrix(as.numeric(onco_result_fixed_5p[i,10+n.fixed+(1:n.fixed^2)]),n.fixed,n.fixed)
  if(det(infor)==0){
    p_sub1[temp] = 1
    
    p_sub2[temp] = 1
  }else{
    score.fix = t(score[1:n.support])
    infor.fix = infor[1:n.support,1:n.support]
    score.random = t(as.numeric(onco_result_casecase_5p[i,10+(1:n.test)]))
    infor.random = t(matrix(as.numeric(onco_result_casecase_5p[i,10+(n.test)+(1:n.test^2)]),
                            n.test,n.test))
    
    p_sub1[temp] = DisplayFixedScoreTestResult(score,infor)
    
    p_sub2[temp] = DisplayMixedScoreTestResult(score.fix,infor.fix,score.random,infor.random)[1]
  }
  
  temp = temp + 1
}
p_sub = list(p_sub1,p_sub2)
save(p_sub,file=paste0("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/p_sub_two_interaction",i1,".Rdata"))
