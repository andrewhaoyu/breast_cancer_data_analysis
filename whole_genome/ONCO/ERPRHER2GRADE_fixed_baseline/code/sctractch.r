#setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result')

i1.idx = c(352,349,347,476,532,209,210,234,247,244)
second.num <- 5
num.of.tumor <- 4

for(j in 1:length(i1.idx)){
  
  i1 <- i1.idx[j]
  print(i1)
  setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result')
  num.total <- 0
  for(i in 1:1000){
    if(i%%100==0){
      print(i)  
    }
   
    load(paste0('ERPRHER2Grade_fixed_onco_resubmit.Rdata',i1,'_',i))
    temp <- length(result.list[[1]])
    num.total <- num.total + temp
  }
  # result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all,score_result_baseline = score_result_baseline,
  #                infor_result_baseline=infor_result_baseline)
  
  score_result <- matrix(0.1,num.total,num.of.tumor+1)
  infor_result <- matrix(0.1,(num.of.tumor+1)*num.total,num.of.tumor+1)
  snpid_result <- rep("c",num.total)
  score_result_baseline <- rep(0,num.total)
  infor_result_baseline <- rep(0,num.total)
  freq.all <- rep(0,num.total)
  num.total <- 0
  for(i in 1:1000){
    if(i%%100==0){
      print(i)  
    }
    
    #print(i)
    load(paste0('ERPRHER2Grade_fixed_onco_resubmit.Rdata',i1,'_',i))
    temp <- length(result.list[[1]])
    snpid_result[num.total+(1:temp)] <- result.list[[1]]
    score_result[num.total+(1:temp),] <- result.list[[2]]
    infor_result[second.num*num.total+(1:(temp*second.num)),] <- result.list[[3]]
    freq.all[num.total+(1:temp)] <- result.list[[4]]
    score_result_baseline[num.total+(1:temp)] <- result.list[[5]]
    infor_result_baseline[num.total+(1:temp)] <- result.list[[6]]
    
    num.total <- num.total + temp
    
  }
  
  result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all,score_result_baseline = score_result_baseline,
                 infor_result_baseline=infor_result_baseline)
  setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
  save(result,file=paste0("./whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2Grade_fixed_onco",i1))
  
}






setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result')
i1 = 58
second.num <- 5
num.of.tumor <- 4
num.total <- 0
for(i in 1:1000){
  print(i)
  load(paste0('ERPRHER2Grade_fixed_baseline_resubmit',i1,'_',i))
  temp <- length(result.list[[1]])
  num.total <- num.total + temp
}

score_result <- matrix(0.1,num.total,num.of.tumor+1)
infor_result <- matrix(0.1,(num.of.tumor+1)*num.total,num.of.tumor+1)
snpid_result <- rep("c",num.total)
score_result_baseline <- rep(0,num.total)
infor_result_baseline <- rep(0,num.total)
freq.all <- rep(0,num.total)
num.total <- 0
for(i in 1:1000){
  print(i)
  load(paste0('ERPRHER2Grade_fixed_baseline_resubmit',i1,'_',i))
  temp <- length(result.list[[1]])
  snpid_result[num.total+(1:temp)] <- result.list[[1]]
  score_result[num.total+(1:temp),] <- result.list[[2]]
  infor_result[second.num*num.total+(1:(temp*second.num)),] <- result.list[[3]]
  freq.all[num.total+(1:temp)] <- result.list[[4]]
  score_result_baseline[num.total+(1:temp)] <- result.list[[5]]
  infor_result_baseline[num.total+(1:temp)] <- result.list[[6]]
  
  num.total <- num.total + temp
  
}

result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all,score_result_baseline = score_result_baseline,
               infor_result_baseline=infor_result_baseline)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
save(result,file=paste0("./whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2Grade_fixed_baseline",i1))


setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result')
i1 = 59
second.num <- 5
num.of.tumor <- 4
num.total <- 0
for(i in 1:1000){
  print(i)
  load(paste0('ERPRHER2Grade_fixed_baseline_resubmit',i1,'_',i))
  temp <- length(result.list[[1]])
  num.total <- num.total + temp
}

score_result <- matrix(0.1,num.total,num.of.tumor+1)
infor_result <- matrix(0.1,(num.of.tumor+1)*num.total,num.of.tumor+1)
snpid_result <- rep("c",num.total)
score_result_baseline <- rep(0,num.total)
infor_result_baseline <- rep(0,num.total)
freq.all <- rep(0,num.total)
num.total <- 0
for(i in 1:1000){
  print(i)
  load(paste0('ERPRHER2Grade_fixed_baseline_resubmit',i1,'_',i))
  temp <- length(result.list[[1]])
  snpid_result[num.total+(1:temp)] <- result.list[[1]]
  score_result[num.total+(1:temp),] <- result.list[[2]]
  infor_result[second.num*num.total+(1:(temp*second.num)),] <- result.list[[3]]
  freq.all[num.total+(1:temp)] <- result.list[[4]]
  score_result_baseline[num.total+(1:temp)] <- result.list[[5]]
  infor_result_baseline[num.total+(1:temp)] <- result.list[[6]]
  
  num.total <- num.total + temp
  
}

result <- list(snpid_reuslt=snpid_result,score_result=score_result,infor_result=infor_result,freq.all=freq.all,score_result_baseline = score_result_baseline,
               infor_result_baseline=infor_result_baseline)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
save(result,file=paste0("./whole_genome/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2Grade_fixed_baseline",i1))




