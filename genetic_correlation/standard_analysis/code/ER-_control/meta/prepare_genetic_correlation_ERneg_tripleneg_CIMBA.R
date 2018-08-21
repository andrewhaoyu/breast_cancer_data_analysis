#################generate the ER- data for genetic correlation analysis
load(paste0("/data/zhangh24/breast_cancer/standard_gwas/ER-_control/meta/ER-_control_meta_result_final.Rdata"))
head(ERNeg_control_meta_result_final)

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/meta_result_shared_1p.Rdata")
head(meta_result_shared_1p)



ERNeg_control_meta_result_final_shared <- merge(ERNeg_control_meta_result_final,
                                                meta_result_shared_1p,
                                                by.x = "rs_id",
                                                by.y = "rs_id")
idx <- which(ERNeg_control_meta_result_final_shared$info>=0.3)
ERNeg_control_meta_result_final_shared <- ERNeg_control_meta_result_final_shared[idx,]
##################

alleles1 <- rep("c",total)
alleles2 <- rep("c",total)
alleles.split.icog <- strsplit(alleles.ICOG,split=":")

alleles.ONCO <- as.character(ICOG.result.clean$SNP.ONCO)
alleles3 <- rep("c",total)
alleles4 <- rep("c",total)
alleles.split.onco <- strsplit(alleles.ONCO,split=":")


for(i in 1:total){
  print(i)
  alleles1[i] <- alleles.split.icog[[i]][3]
  alleles2[i] <- alleles.split.icog[[i]][4]
  alleles3[i] <- alleles.split.onco[[i]][3]
  alleles4[i] <- alleles.split.onco[[i]][4]
}

allele_temp <- ERNeg_control_meta_result_final_shared$var_name
ERNeg_control_meta_result <- ERNeg_control_meta_result_final_shared[,c(1,2,3,4,5,)]