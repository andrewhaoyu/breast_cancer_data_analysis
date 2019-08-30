# args = commandArgs(trailingOnly = T)
# i1 = as.numeric(args[[1]])
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared_1p_sex.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared_1p_sex.Rdata")
#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared.Rdata")
#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared.Rdata")
# load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_only_shared_1p_082119.Rdata")
# load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_only_shared_1p_082119.Rdata")

second.num <- 5

icog_result_shared_1p_casecase <- icog_result_shared_1p
onco_result_shared_1p_casecase <- onco_result_shared_1p
# icog_result_only_shared_1p_casecase <- icog_result_only_shared_1p
# onco_result_only_shared_1p_casecase <- onco_result_only_shared_1p


icog_score_infor_casecase <- icog_result_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]
onco_score_infor_casecase <- onco_result_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]
# 
# icog_score_infor_icog_only_casecase <- icog_result_only_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]
# 
# onco_score_infor_icog_only_sudo_casecase <- icog_score_infor_icog_only_casecase
# onco_score_infor_icog_only_sudo_casecase[] <- 0
# for(i in 1:nrow(onco_score_infor_icog_only_sudo_casecase)){
#   onco_score_infor_icog_only_sudo_casecase[i,(second.num+1):(second.num+second.num^2)] <- as.vector(diag(10000,second.num))
# }
# onco_score_infor_onco_only_casecase <- onco_result_only_shared_1p_casecase[,11:(11+second.num+second.num^2-1)]
# icog_score_infor_onco_only_sudo_casecase <- onco_score_infor_onco_only_casecase
# icog_score_infor_onco_only_sudo_casecase[] <- 0
# for(i in 1:nrow(icog_score_infor_onco_only_sudo_casecase)){
#   icog_score_infor_onco_only_sudo_casecase[i,(second.num+1):(second.num+second.num^2)] <- as.vector(diag(10000,second.num))
# }
# 

icog_onco_score_infor_casecase <- cbind(icog_score_infor_casecase,onco_score_infor_casecase)
# icog_onco_score_infor_icog_only_casecase <- cbind(icog_score_infor_icog_only_casecase,onco_score_infor_icog_only_sudo_casecase)
# icog_onco_score_infor_onco_only_casecase <- cbind(icog_score_infor_onco_only_sudo_casecase,onco_score_infor_onco_only_casecase)

# icog_onco_score_infor_final_casecase <- rbind(icog_onco_score_infor_casecase,icog_onco_score_infor_icog_only_casecase,icog_onco_score_infor_onco_only_casecase)
# icog_onco_score_infor_casecase <- icog_onco_score_infor_final_casecase

n <- nrow(icog_onco_score_infor_casecase)
# debug.one.line <- c(rep(0,second.num),as.vector(diag(second.num)),rep(0,second.num),as.vector(diag(second.num)))
# debug.idx <- c(9649548,9650051)
# for(k in 1:length(debug.idx)){
#   print(k)
#   icog_onco_score_infor_casecase[debug.idx[k],] <- debug.one.line  
# }
# 

icog_onco_score_infor <- icog_onco_score_infor_casecase



library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")

rm(icog_result_shared_1p)
rm(onco_result_shared_1p)
rm(icog_result_shared_1p_casecase)
rm(onco_result_shared_1p_casecase)
gc()

#library(foreach)
#library(doParallel)
#no.cores <- 12
n <- nrow(icog_onco_score_infor)
#icog_onco_score_infor_temp <- icog_onco_score_infor[1:10^5,]
#n <- nrow(icog_onco_score_infor_temp)
#pvalue <- rep(0,n)
#size <- 1000
#n <- nrow(icog_onco_score_infor)


# fixed.second.num <- 2
# random.second.num <- 3
#registerDoParallel(no.cores)
#pvalue <- rep(0,nrow(icog_onco_score_infor))


#pvalue.list <- foreach(i=1:size)%dopar%
#{
#print(i)
# start.end <- startend(n,size,i1)
# start <- start.end[1]
# end <- start.end[2]
start <- 1
end <- nrow(icog_onco_score_infor)
result_summary <- matrix(0,end-start+1,(second.num+second.num^2+1))
#pvalue_sub <- rep(0,end-start+1)
temp = 1
for(j in start:end){
  if(j%%1000==0){
    print(j)  
  }
  
  icog_onco_score_infor_oneline <- icog_onco_score_infor[j,]
  #some SNPs were put infor and score are 0 because of unconvergence
  icog_infor <- matrix(icog_onco_score_infor_oneline[(second.num+1):(second.num+second.num^2)],second.num,second.num)
  onco_infor <- matrix(icog_onco_score_infor_oneline[(second.num+second.num^2+second.num+1):(second.num+second.num^2+second.num+second.num^2)],second.num,second.num)
  if(det(icog_infor)==0){
    print(j)
    icog_onco_score_infor_oneline[(second.num+1):(second.num+second.num^2)] <- as.vector(diag(10000,second.num))
  }
  if(det(onco_infor)==0){
    print(j)
    icog_onco_score_infor_oneline[(second.num+second.num^2+second.num+1):(second.num+second.num^2+second.num+second.num^2)] <- as.vector(diag(10000,second.num))
  }
  
  result_temp <- MetaFixedPfunction_temp(icog_onco_score_infor_oneline,second.num)
  result_summary[temp,1:second.num] <- as.vector(result_temp[[1]])
  result_summary[temp,(second.num+1):(second.num+second.num^2)] <- as.vector(result_temp[[2]])
  result_summary[temp,(second.num+second.num^2+1)] <- as.numeric(result_temp[[3]][11])
  temp = temp+1
}

meta_result_shared_1p <- icog_result_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]


meta_result_shared_1p_no_pvalue <- rbind(meta_result_shared_1p,meta_result_shared_1p_icog_only,meta_result_shared_1p_onco_only)


p.value <- result_summary[,31]
log.odds <- result_summary[,1:5]
sigma <- result_summary[,6:30]


meta_result_shared_1p <- cbind(meta_result_shared_1p_no_pvalue,p.value,log.odds,sigma)

save(meta_result_shared_1p,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_sex.Rdata"))



colnames(meta_result_shared_1p)[21:45] <- 
  paste0("cov",c(1:25))



library(data.table)
#load CIMBA BCAC meta-analysis data
#chromosome 23 data are not meta-analyzed yet
CIMBA.result <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/brca1_bcac_tn_meta.txt"))



CIMBA.BCAC <- merge(meta_result_shared_1p,CIMBA.result,by.x = "var_name",by.y = "MarkerName")
idx.am <- which(CIMBA.BCAC$a1!=toupper(CIMBA.BCAC$Allele2))
#reverse the allele that not aligned in bcac and CIMBA
temp1 = CIMBA.BCAC$a1[idx.am] 
temp2 = CIMBA.BCAC$a0[idx.am]
CIMBA.BCAC$a1[idx.am] = temp2
CIMBA.BCAC$a0[idx.am] = temp1

CIMBA.BCAC$Effect[idx.am] <- -CIMBA.BCAC$Effect[idx.am]
CIMBA.BCAC$Freq1[idx.am] <- 1-CIMBA.BCAC$Freq1[idx.am]
CIMBA.BCAC$MinFreq[idx.am] <- 1-CIMBA.BCAC$MinFreq[idx.am]
CIMBA.BCAC$MaxFreq[idx.am] <- 1-CIMBA.BCAC$MaxFreq[idx.am]

CIMBA.BCAC_update <- CIMBA.BCAC %>% 
  mutate(var_effect = StdErr^2) %>% 
  select(var_name,"5",cov25,Effect,var_effect)

#BCAC TN and CIMBA BRCA1 meta analysis
meta_effect <- rep(0,nrow(CIMBA.BCAC_update))
meta_std <- rep(0,nrow(CIMBA.BCAC_update))
meta_p <- rep(0,nrow(CIMBA.BCAC_update))
for(i in 1:nrow(CIMBA.BCAC_update)){
  if(i%%1000==0){
    print(i)  
  }
  
  result_temp <- MetaFixedPfunction_temp(CIMBA.BCAC_update[i,2:5],1)
  meta_effect[i] <- result_temp[[1]]
  meta_std[i] <- sqrt(result_temp[[2]])
  meta_p[i] <- as.numeric(result_temp[[3]][2])
  
}



CIMBA.BCAC_update$meta_effect <- meta_effect
CIMBA.BCAC_update$meta_std <- meta_std
CIMBA.BCAC_update$meta_p <- meta_p
#Some alleles change the order due to BCAC and CIMBA meta-analysis.
#we need to change them back
colnames(CIMBA.BCAC_update)[1] <- "MarkerName"
CIMBA.BCAC_update$a1_update <- CIMBA.BCAC$a1
CIMBA.BCAC_update$a0_update <- CIMBA.BCAC$a0
CIMBA.BCAC_update$Freq1_update <- CIMBA.BCAC$Freq1
CIMBA.BCAC_update$MinFreq_update <- CIMBA.BCAC$MinFreq
CIMBA.BCAC_update$MaxFreq_update <- CIMBA.BCAC$MaxFreq





#update CIMBA sex chromosome data
CIMBA.result.update <- left_join(CIMBA.result,CIMBA.BCAC_update,by ="MarkerName")
idx <- which(!is.na(CIMBA.result.update$Effect.y))
#use the update information allele1,allele2,freq1,minfreq,maxfreq,effect.x, stderr,p-value
CIMBA.result.update[idx,c(2,3,4,6,7,8,9,10)] <-
  CIMBA.result.update[idx,c(26,25,28,27,29,20,23,24)]
head(CIMBA.result.update)
#drop the uncessary data
CIMBA.result.update <- CIMBA.result.update[,c(1:17)]
write.table(CIMBA.result.update,file = "/spin1/users/zhangh24/breast_cancer_data_analysis/data/CIMBA_BCAC_meta_analysis_083019.txt",quote=F,row.names=F)









