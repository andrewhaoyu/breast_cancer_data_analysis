#-------------------------------------------------------------------
# Update Date: 11/19/2018
# Create Date: 11/19/2018
# Goal: Take out two-stage model results for HER2 
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#Load mixed effect model results
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata")
mixed <- meta_result_shared_1p
#Load fixed effect model results
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
fixed <- meta_result_shared_1p
#Load intirnsic subtypes results
load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p.Rdata"))
intrin <- meta_result_shared_1p
# all.equal(mixed$rs_id,fixed$rs_id)
# all.equal(mixed$rs_id,intrin$rs_id)
idx <- which(mixed$CHR==11&
               mixed$position>=74459913
             &mixed$position<= 74553458)
snp.info <- mixed[,c(1:14)]

new.data <- 