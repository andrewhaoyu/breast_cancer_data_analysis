#meta-analysis for icogs and oncoarray result
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
library(Rcpp)
library(RcppArmadillo)
sourceCpp("/gpfs/gsfs11/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/code/ICOG/Saddle.cpp")

#MetaPfunction can perform the meta-analysis for icogs and oncoarray 
MetaPfunction <- function(icog_onco_score_infor_one, second.num){
  #get the score for icogs
  score.icog <- as.numeric(icog_onco_score_infor_one[1:(second.num)])
  #get the information matrix for icogs
  infor.icog <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1):
                                                              (second.num+second.num^2) ]),
                       ncol = second.num)
  start <- second.num+second.num^2
  #get the score for onco
  score.onco <- as.numeric(icog_onco_score_infor_one[(1+start):
                                                       (second.num+start)])
  #get the information matrix for onco
  infor.onco <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1+start):
                                                              (second.num+second.num^2+start) ]),ncol=second.num)
  
  
  meta.result <- ScoreMetaAnalysis(score.icog,infor.icog,
                                   score.onco,infor.onco)
  score.meta <- t(meta.result[[1]])
  infor.meta <- meta.result[[2]]
  #calculate the p-value for based on score statistics and information matrix
  #to perform fixed-effect test using chi-square test
  #we can use function ScoreGlobalTestForAssoc(score, infor)
  #to perform random-effect test using mixture chi-square test
  #we can use function ScoreMixedGlobalTestForHeter(score, infor)
  p.value = ScoreGlobalTestForAssoc(score.meta, infor.meta)
  return(p.value)
}

ScoreGlobalTestForAssoc <- function(score,infor){
  infor <- as.matrix(infor)
  df <- length(score)
  GTA.stat <- score%*%solve(infor)%*%t(score)
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  places = 3
  power.number <- floor(-log10(p.value.GTA))+places
  if(p.value.GTA==0){
    return(1E-300)
  }else{
    ###format the output with three digits in total
    p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)
    
  }
  
  return(p.value.GTA)
  
}


#score statistics meta-analysis
ScoreMetaAnalysis <- function(score1,infor1,score2,infor2){
  infor.meta <- infor1+infor2
  score.meta <- score1+score2
  
  return(list(score.meta = score.meta,
              infor.meta = infor.meta))
  
}

#Random effect model score test
#Q_simga = t(score)%*%score
#test statistics follow a mixture of chi-square test statistics
#the weights of the mixture chi-square statistics are the eigen values of information matrix
ScoreMixedGlobalTestForHeter <- function(score.casecase,infor.casecase){
  
  
  GTH.stat <- as.numeric(score.casecase%*%t(score.casecase))
  lamda <- eigen(infor.casecase)$values
  
  acc = 1e-09
  lim = 2000000
  
  #davies method is a numerical algorithm
  #lim is the number of monte carol runs
  #acc is the accuracy
  #you can increase the lim number and decrease the acc level if the function speed is really fast
  result <- davies(GTH.stat,lamda,lim = lim,acc=acc)
  p.value.GTH <- result[[3]]
  
  if(result[[2]]!=0){
    #if result[[2]] is not 0
    #davies method didn't converge. 
    #I noticed this problem happens for 6_102485159_G_A and 1_87277974_G_A;
    #These two SNPs are oncoarray only SNPs with allele frequency around 0.01
    #just put the p-value as 1 for these SNPs
    print("chisq p value accuracy could't be reached")
    p.value.GTH = 1
  }
  
  if(p.value.GTH <0){
    #davies method is a numerical algorithm
    #when p-value is really small, the value can be negative
    #if it happens, just put the p-value as 1e-09 which is the accuracy for the function
    p.value.GTH <- acc
  }
  #places <- 3
  #power.number <- floor(-log10(p.value.GTH))+places
  #p.value.GTH <- round(p.value.GTH*10^power.number)/(10^power.number)
  
  
  return(p.value.GTH)
  
  
  
}


ScoreMixedGlobalTestForHeterUpdate <- function(score.casecase,infor.casecase){
  
  
  GTH.stat <- as.numeric(score.casecase%*%t(score.casecase)) #score*t(score)
  lambda <- eigen(infor.casecase)$values #the eigen value of information matrix
  
  p.value.GTH <- Saddle(GTH.stat,lambda) #p-value for linear combination of chi-square
  
  if(p.value.GTH==2){
    acc = 1e-09
    lim = 2000000
    
    #davies method is a numerical algorithm
    #lim is the number of monte carol runs
    #acc is the accuracy
    #you can increase the lim number and decrease the acc level if the function speed is really fast
    result <- davies(GTH.stat,lambda,lim = lim,acc=acc)
    p.value.GTH <- result[[3]]
    
    if(result[[2]]!=0){
      #if result[[2]] is not 0
      #davies method didn't converge. 
      #I noticed this problem happens for 6_102485159_G_A and 1_87277974_G_A;
      #These two SNPs are oncoarray only SNPs with allele frequency around 0.01
      #just put the p-value as 1 for these SNPs
      print("chisq p value accuracy could't be reached")
      p.value.GTH = 1
    }
    
    if(p.value.GTH <0){
      #davies method is a numerical algorithm
      #when p-value is really small, the value can be negative
      #if it happens, just put the p-value as 1e-09 which is the accuracy for the function
      p.value.GTH <- acc
    }
    #places <- 3
    #power.number <- floor(-log10(p.value.GTH))+places
    #p.value.GTH <- round(p.value.GTH*10^power.number)/(10^power.number)
    
    
   
  }
  
  return(p.value.GTH)
  
  
  
}
#install R package
#bc2 is a development version of TOP package
#I used bc2 in my previous analyses
#the function of bc2 and TOP are almost the same
#TOP has more documentation
#to install bc2 or TOP, one needs to use install_github function
#you can specify the directory to your local directory
#you need bc2 package to run the analysis in this code
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/4.2/", install_github('andrewhaoyu/bc2'))
library(bc2, 
        lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/4.2/")
library(CompQuadForm)
#specify the second.num given the model
#additive model second. num is 5
second.num <- 5
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/icog_onco_score_infor.rdata")
#count the total number of SNPs
n <- nrow(icog_onco_score_infor)
icog_onco_score_infor[1028143,]

# #SNPs 6_102485159_G_A and 1_87277974_G_A are oncoarray only SNPs with allele frequency around 0.01
# #during numerical 
# debug.one.line <- c(rep(0,second.num),as.vector(diag(second.num)),rep(0,second.num),as.vector(diag(second.num)))
# debug.idx <- c(9649548,9650051)
# icog_onco_score_infor_final[9649548,]
# for(k in 1:length(debug.idx)){
# print(k)
#     icog_onco_score_infor[debug.idx[k],] <- debug.one.line  
# }

#split the jobs into 1000 subjobs
size <- 1000
#run the job from start to end
start.end <- startend(n,size,i1)
start <- start.end[1]
end <- start.end[2]
#count the p_value
pvalue_sub <- rep(0,end-start+1)
temp = 1

for(j in start:end){
  print(j)
  
  icog_onco_score_infor_oneline <- icog_onco_score_infor[j,]
  pvalue_sub[temp] <- MetaPfunction(icog_onco_score_infor_oneline,second.num)
  temp = temp+1
}


#save the pvalue_sub file to a folder
save(pvalue_sub,file=paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/p_value_sub",i1,".Rdata"))
# 
# idx = c(1028143, 1338930, 1579924, 2290249, 2290288, 2290311, 3211390, 3796181, 4690893, 5041247, 5041251, 5077860, 5077863, 5077864, 5077865, 5173984, 5360746, 5935978, 5935983, 5935985, 5936002, 5936005, 5936007, 5936008, 5936016, 5936023, 6049450, 6049451, 7677446)
# pvalue_sub = length(idx)
# temp = 1
# for(j in idx){
#   print(j)
#   
#   icog_onco_score_infor_oneline <- icog_onco_score_infor[j,]
#   pvalue_sub[temp] <- MetaPfunction(icog_onco_score_infor_oneline,second.num)
#   temp = temp+1
# }
# load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_final.Rdata")
# meta_result_shared_1p[idx,"p.value"] = 1E-300
# save(meta_result_shared_1p, file = "/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_final.Rdata")
#meta_result_shared_1p[idx,1:14]
idx <- which(data$var_name == "2_172376262_A_C")
data[idx,]
