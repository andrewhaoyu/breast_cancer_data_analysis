load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared_1p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared_1p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_only_shared_1p.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_only_shared_1p.Rdata")


meta_result_shared_1p <- icog_result_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]
meta_result_shared_1p_icog_only <- icog_result_only_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]
meta_result_shared_1p_onco_only <- onco_result_only_shared_1p[,c(1:10,(ncol(icog_result_shared_1p)-3):ncol(icog_result_shared_1p))]


meta_result_shared_1p_no_pvalue <- rbind(meta_result_shared_1p,meta_result_shared_1p_icog_only,meta_result_shared_1p_onco_only)


# save(meta_result_shared_1p_no_pvalue,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_no_pvalue.Rdata")


n <- nrow(meta_result_shared_1p_no_pvalue)

p.value <- rep(0,n)
total <- 0
log.odds <- matrix(0,n,5)
sigma <- matrix(0,n,25)
for(i1 in 1:1000){
print(i1)  
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/reuslt_summary_sub",i1,".Rdata"))
  temp <- nrow(result_summary)
  
  p.value[total+(1:temp)] <- result_summary[,31]
  log.odds[total+(1:temp),] <- result_summary[,1:5]
  sigma[total+(1:temp),1:25] <- result_summary[,6:30]
  total <- total + temp
  
}


meta_result_shared_1p <- cbind(meta_result_shared_1p_no_pvalue,p.value,log.odds,sigma)
save(meta_result_shared_1p,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p.Rdata"))

fine_mapping <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/fine_mapping_regions.csv",header= T)

idx <- which(intrinsic_subtype_triple_negative_results$position==7571752&intrinsic_subtype_triple_negative_results$CHR==17)
intrinsic_subtype_triple_negative_results[idx,]

intrinsic_subtype_triple_negative_results <- meta_result_shared_1p[,c(2,13,14,3,11,4,5,20,45)]
colnames(intrinsic_subtype_triple_negative_results) <- 
  c("rs_id","SNP_ICOGs","SNP_ONCO","position",
    "CHR","freq_a1","imputation_quality","log_odds_triple_negative","variance_triple_negative")
save(intrinsic_subtype_triple_negative_results,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/intrinsic_subtype_triple_negative_results.Rdata"))
alleles.ICOG <- as.character(intrinsic_subtype_triple_negative_results$SNP_ICOGs)
total <- nrow(intrinsic_subtype_triple_negative_results)
alleles1 <- rep("c",total)
alleles2 <- rep("c",total)
alleles.split.icog <- strsplit(alleles.ICOG,split=":")

alleles.ONCO <- as.character(intrinsic_subtype_triple_negative_results$SNP_ONCO)
alleles3 <- rep("c",total)
alleles4 <- rep("c",total)
alleles.split.onco <- strsplit(alleles.ONCO,split=":")


for(i in 1:total){
  if(i%%100==0){
    print(i)  
  }
  alleles1[i] <- alleles.split.icog[[i]][3]
  alleles2[i] <- alleles.split.icog[[i]][4]
  alleles3[i] <- alleles.split.onco[[i]][3]
  alleles4[i] <- alleles.split.onco[[i]][4]
}

alleles.data <- data.frame(alleles1,alleles2,alleles3,alleles4)





idx <- which(is.na(alleles1)&!is.na(alleles3))
alleles1[idx] <- alleles3[idx]
alleles2[idx] <- alleles4[idx]
head(intrinsic_subtype_triple_negative_results)


idx <- which(is.na(alleles1))
head(alleles.data[idx,])

alleles.ONCO <- as.character(intrinsic_subtype_triple_negative_results$SNP_ONCO[idx])
alleles.split.onco <- strsplit(alleles.ONCO,split="_")
for(i in 1:length(idx)){
  # if(i%%100==0){
  #   print(i)  
  # }
  alleles1[idx[i]] <- alleles.split.onco[[i]][3]
  alleles2[idx[i]] <- alleles.split.onco[[i]][4]
  alleles3[idx[i]] <- alleles.split.onco[[i]][3]
  alleles4[idx[i]] <- alleles.split.onco[[i]][4]
}
alleles.data <- data.frame(alleles1,alleles2,alleles3,alleles4)
intrinsic_subtype_triple_negative_results$allele1 <- alleles.data$alleles1
intrinsic_subtype_triple_negative_results$allele2 <- alleles.data$alleles2
idx <- which(is.na(alleles1))
length(idx)
head(intrinsic_subtype_triple_negative_results[idx,])
head(alleles.data[idx,])
intrinsic_subtype_triple_negative_results[idx[1:100],]
save(intrinsic_subtype_triple_negative_results,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/intrinsic_subtype_triple_negative_results.Rdata"))
# (intrinsic_subtype_triple_negative_results[idx[1:10],])



idx_cut <- NULL
start <- fine_mapping$start
end <- fine_mapping$end
CHR <- fine_mapping$V3


for(i in 1:nrow(fine_mapping)){
  print(i)
  chr_temp <- CHR[i]
  start_temp <- start[i]
  end_temp <- end[i]
  idx <- which(meta_result_shared_1p$CHR==chr_temp&meta_result_shared_1p$position>=start_temp&
                 meta_result_shared_1p$position<=end_temp)
  idx_cut <- c(idx_cut,idx)
}
############duplicate variables won't mater
idx_cut <- unique(idx_cut)
meta_result_shared_1p_filter <- meta_result_shared_1p[-idx_cut,]

save(meta_result_shared_1p_filter,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_filter_1M.Rdata")


new_filter <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
new_filter <- new_filter[1:22,]

new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))
new_filter2 <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Filter_based_on_two_stage.csv",header=T,stringsAsFactors = F)
new_filter2 <- new_filter2[1:19,1:3]

new_filter <- rbind(new_filter,new_filter2)
idx_cut <- NULL

position.cut <- 500*10^3

for(i in 1:nrow(new_filter)){
  print(i)
  chr_temp <- new_filter[i,3]
  position_temp <- new_filter[i,2]
  position_low <- position_temp-position.cut
  position_high <- position_temp+position.cut
  idx <- which(meta_result_shared_1p_filter$CHR==chr_temp&meta_result_shared_1p_filter$position>position_low&
                 meta_result_shared_1p_filter$position<position_high)
  idx_cut <- c(idx_cut,idx)
}
idx_cut <- unique(idx_cut)
meta_result_shared_1p_filter_Ju <- meta_result_shared_1p_filter[-idx_cut,]

idx.min <- which(meta_result_shared_1p_filter_Ju$p.value<=5E-08)
sig_result <- meta_result_shared_1p_filter_Ju[idx.min,]

idx.test <- which((sig_result$position%in%extract.list.ld$position)&(sig_result$CHR%in%extract.list.ld$CHR))


save(meta_result_shared_1p_filter_Ju,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_filter_1M_Ju.Rdata")



load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata"))



meta_result_shared_1p <- meta_result_shared_1p[,c(2,11,3,15)]
write.table(meta_result_shared_1p,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_mixed.txt"),row.names = F,quote = F)







# idx <- which.min(meta_result_shared_1p_filter_Ju$p.value)
# 
# meta_result_shared_1p_filter_Ju[idx,]





# load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M_Ju.Rdata")
# 
# idx <- which(meta_result_shared_1p_filter_Ju$CHR==1&meta_result_shared_1p_filter_Ju$position== 120485335)
# meta_result_shared_1p_filter_Ju[idx,]
# 
# idx <- which(meta_result_shared_1p_filter_Ju$CHR==2&meta_result_shared_1p_filter_Ju$position== 67902524)
# meta_result_shared_1p_filter_Ju[idx,]
# idx <- which(meta_result_shared_1p_filter_Ju$CHR==11&meta_result_shared_1p_filter_Ju$position== 120233626)
# meta_result_shared_1p_filter_Ju[idx,]
# idx <- which(meta_result_shared_1p_filter_Ju$CHR==18&meta_result_shared_1p_filter_Ju$position== 10354649)
# meta_result_shared_1p_filter_Ju[idx,]
# 
# 


























