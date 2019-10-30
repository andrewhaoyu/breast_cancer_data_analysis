setwd("/data/zhangh24/breast_cancer_data_analysis/")
num.total <- 0
for(i1 in 232){
  for(i2 in 1:1000){
    print(i2)
    load(paste0("./whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/ERPRHER2Grade_casecase_onco_resubmit",i1,"_",i2))
    temp <- length(result[[1]])
    num.total <- temp+num.total
  }
}

snpid <- rep("c",num.total)
score <- matrix(0,num.total,3)
infor <- matrix(0,num.total,9)
freq <- rep(0,num.total)
num.total <- 0
for(i1 in 232){
  for(i2 in 1:1000){
    print(i2)
    load(paste0("./whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/ERPRHER2Grade_casecase_onco_resubmit",i1,"_",i2))
    temp <- length(result[[1]])
    snpid[num.total+(1:temp)] <- result[[1]]
    score[num.total+(1:temp),] <- result[[2]]
    infor[num.total+(1:temp),] <- result[[3]]
    freq[num.total+(1:temp)] <- result[[4]]
    num.total <- temp+num.total
  }
}

result <- list(snpid,score,infor,freq)
save(result,file=paste0("./whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/ERPRHER2Grade_casecase_onco",i1))
