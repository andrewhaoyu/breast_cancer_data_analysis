setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
load(paste0("./discovery_SNP/CIMBA_conditional_analysis/result/CIMBA_snp_name_match.Rdata"))
load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/intrinsic_subtype_triple_negative_results.Rdata"))
log.odds <- rep(0,4121)
var.log.odds <- rep(0,4121)
total <- 0
for(i1 in 1:825){
  print(i1)
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/conditional_result",i1,".Rdata"))
  temp <- nrow(conditional_result)
  log.odds[total+(1:temp)] <- conditional_result[,5]
  var.log.odds[total+(1:temp)] <- conditional_result[,30]
  total <- total+temp
}
chr.pos <- paste0(CIMBA_snp$CHR,":",CIMBA_snp$position)
conditional_result <- data.frame(log.odds,var.log.odds,CIMBA_snp$rs_id)
colnames(conditional_result)[3] <- "rs_id"
head(intrinsic_subtype_triple_negative_results)
head(conditional_result)
# intrinsic_subtype_triple_negative_results$chr.pos <- paste0(intrinsic_subtype_triple_negative_results$CHR,":",
#                                                             intrinsic_subtype_triple_negative_results$position)

new.data <- merge(conditional_result,intrinsic_subtype_triple_negative_results,by.x = "rs_id",by.y = "rs_id")

intrinsic_subtype_triple_negative_results_chr12 <- 
  new.data[,c(1,4,5,6,7,8,9,2,3,12,13)]
colnames(intrinsic_subtype_triple_negative_results_chr12)[c(8,9)] <- c("log_odds_triple_negative","variance_triple_negative")
save(intrinsic_subtype_triple_negative_results_chr12,
     file = "/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/CIMBA_conditional_analysis/result/intrinsic_subtype_triple_negative_results_chr12.Rdata")






# try <- match(conditional_result$chr.pos,intrinsic_subtype_triple_negative_results$chr.pos)
# 
# try.data <- intrinsic_subtype_triple_negative_results[try,]
# 
# 
# try2 <- try.data$chr.pos%in%new.data$chr.pos
