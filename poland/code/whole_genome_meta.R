arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
load("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_fixed_5p.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/onco_result_casecase_5p.Rdata")
library(bc2)
num = nrow(onco_result_fixed_5p)
size = 1000
start.end <- startend(num,size,i1)
start <- start.end[1]
end <- start.end[2]
file.num <- end-start+1
p_sub1 <- rep(0,file.num)
p_sub2 <- rep(0,file.num)
temp = 1
for(i in start:end){
    score = t(as.numeric(onco_result_fixed_5p[i,11:15]))
  infor = matrix(as.numeric(onco_result_fixed_5p[i,16:40]),5,5)
  p_sub1[temp] = DisplayFixedScoreTestResult(score,infor)
  score.fix = t(score[1:2])
  infor.fix = infor[1:2,1:2]
  score.random = t(as.numeric(onco_result_casecase_5p[i,11:13]))
  infor.random = t(matrix(as.numeric(onco_result_casecase_5p[i,14:22]),
                        3,3))
  p_sub2[temp] = DisplayMixedScoreTestResult(score.fix,infor.fix,score.random,infor.random)[1]
  temp = temp + 1
}
p_sub = list(p_sub1,p_sub2)
save(p_sub,file=paste0("/data/zhangh24/breast_cancer_data_analysis/poland/result/whole_genome/p_sub",i1,".Rdata"))
