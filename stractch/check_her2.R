#-------------------------------------------------------------------
# Update Date: 11/20/2018
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
snp.info <- mixed[idx,c(1:14)]
mixed.p <- mixed[idx,c(15)]
fixed.p <- fixed[idx,c(15)]
names.subtypes <-  c("Luminal_A","Luminal_B",
                     "Luminal_B_HER2Neg",
                     "HER2Enriched",
                     "TripleNeg")
intrin.log <- intrin[idx,c(16:20)]
intrin.sigma <- intrin[idx,(21:45)]
intrin.result <- NULL
library(bc2)
for(i in 1:nrow(intrin.log)){
   intrin.result <- rbind(intrin.result,
                          DisplaySecondStageTestResult(as.numeric(intrin.log[i,]),
                                                       (matrix(as.numeric(intrin.sigma[i,]),5,5)))) 
}

colnames(intrin.result)[c(1,3,5,7,9)] <- paste0(names.subtypes,"_",c("OR (95%CI)"))
colnames(intrin.result)[c(2,4,6,8,10)] <- paste0(names.subtypes,"_",c("p-value")) 
new.data <- cbind(snp.info,mixed.p,fixed.p,intrin.result[,1:10])
colnames(new.data)[c(15,16)] <- c("Mixed model global association test",
                                  "Fixed model global association test")
write.csv(new.data,file="/spin1/users/zhangh24/breast_cancer_data_analysis/data/HER2_check.csv",quote=F,row.names=F)
