library(data.table)
library(bc2)
#data <- fread("./data/dataset_montse_20180522.txt")
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis')
data <- fread("./data/dataset_montse_20180524.txt")


transform <- function(status,tumor){
  n <- length(status)
  result <- tumor
  idx.mis <- which(status==1&is.na(tumor))
  result[idx.mis] <- 888
  idx.control <- which(status==0)
  result[idx.control] <- NA
  return(result)
}
##########remove status 2 in-situ
##########remove status 3 unknown people
idx.remove <- which(data$status==2|data$status==3)
data <- data[-idx.remove,]
#########remove people with missing parity_cat


idx.remove <- which(data$parity_cat==9)
data <- data[-idx.remove,]


status <- data$status
ER <- data$ER_status1
ER.update <- transform(status,ER)
PR <- data$PR_status1
PR.update <- transform(status,PR)
HER2 <- data$HER2_status1
HER2.update <- transform(status,HER2)
Grade <- data$Grade1
Grade.update <- transform(status,Grade)
table(data$parity_cat)
table(data$study)
table(data$ethnicity)
table(data$refage)
y.pheno.mis <- cbind(status,ER.update,PR.update,HER2.update
                     ,Grade.update)
colnames(y.pheno.mis) <- c("casecon","ER","PR","HER","Grade")
parity.mat <- model.matrix(~as.factor(data$parity_cat)-1)[,-1]
colnames(parity.mat) <- paste0("parity_cat",c(1:4))
study.mat <- model.matrix(~as.factor(data$study)-1)[,-1]
colnames(study.mat) <- unique(data$study)[-1]
ethnicity.mat <- model.matrix(~as.factor(data$ethnicity)-1)[,-1]
colnames(ethnicity.mat) <- paste0("ethnicity",c(1:6))[-1]
refage <- data$refage

# x.covar <- cbind(parity.mat,ethnicity.mat,refage)
# 
# #parity.mat,
# #,
# 
# model <- TwoStageModel(y = y.pheno.mis,
#                        #baselineonly = study.mat,
#                        additive = x.covar,
#                        missingTumorIndicator = 888)
library(xlsx)
# write.xlsx(model[[4]],file="./risk_factor/result/risk_factor_no_study.xlsx",sheetName="second_stage_parameter")
# write.xlsx(model[[5]],file="./risk_factor/result/risk_factor_no_study.xlsx",sheetName="test_result",append=T)
# write.xlsx(model[[7]],file="./risk_factor/result/risk_factor_no_study.xlsx",sheetName="first_stage_result",append=T)

parity <- data$parity_cat
x.covar <- cbind(parity,ethnicity.mat,refage)
model.lin <- TwoStageModel(y = y.pheno.mis,
                       #baselineonly = study.mat,
                       additive = x.covar,
                       missingTumorIndicator = 888)
write.xlsx(model.lin[[4]],file="./risk_factor/result/risk_factor_paraty_lin.xlsx",sheetName="second_stage_parameter")
write.xlsx(model.lin[[5]],file="./risk_factor/result/risk_factor_paraty_lin.xlsx",sheetName="test_result",append=T)
# write.xlsx(model[[7]],file="./risk_factor/result/risk_factor_no_study.xlsx",sheetName="first_stage_result",append=T)


load("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/z.design.Rdata")

model.lin.tninter <- EMmvpolySelfDesign(y.pheno.mis,x.self.design = x.covar[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.covar[,2:ncol(x.covar)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

write.xlsx(model.lin.tninter[[4]],file="./risk_factor/result/risk_factor_paraty_lin.xlsx",sheetName="second_stage_parameter_third_order",append=T)
write.xlsx(model.lin.tninter[[5]],file="./risk_factor/result/risk_factor_paraty_lin.xlsx",sheetName="test_result_third_order",append=T)

z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2-",
                        "HER2 Enriched",
                        "Triple Negative")
design_cat <- data$design_cat
#idx.pop <- which(data$design_cat==0)
parity.mat_design_inter <- design_cat*parity.mat
colnames(parity.mat_design_inter) <- paste0(colnames(parity.mat),"_interaction")
x.covar <- cbind(parity.mat,parity.mat_design_inter,design_cat,ethnicity.mat,refage)
model2 = EMmvpolySelfDesign(y.pheno.mis,x.self.design = x.covar,z.design=z.design,
                            baselineonly = NULL,additive =NULL,pairwise.interaction = NULL,
                            saturated = NULL,
                            missingTumorIndicator = 888)
model2[[4]][,2] <- rep(c("Luminial A","Luminal B",
  "Luminal B HER2-",
  "HER2 Enriched",
  "Triple Negative"),15)
write.xlsx(model2[[4]],file="./risk_factor/result/risk_factor_paraty_lin.xlsx",sheetName="intrinsic_subtye",append=T)

idx.pop <- which(data$design_cat==0)
idx.non_missing <- which(data$molgroup!=9)
data.pop.com <- data[idx.pop[idx.non_missing],]


# parity.mat_design_inter <- design_cat*parity.mat
# colnames(parity.mat_design_inter) <- paste0(colnames(parity.mat),"_interaction")
x.covar <- cbind(parity.mat,ethnicity.mat,refage)
model3 = EMmvpolySelfDesign(y.pheno.mis[idx.pop,],x.self.design = x.covar[idx.pop,],z.design=z.design,
                            baselineonly = NULL,additive =NULL,pairwise.interaction = NULL,
                            saturated = NULL,
                            missingTumorIndicator = 888)
model3[[4]][,2] <- rep(c("Luminial A","Luminal B",
                         "Luminal B HER2-",
                         "HER2 Enriched",
                         "Triple Negative"),10)
write.xlsx(model3[[4]],file="./risk_factor/result/risk_factor_paraty_lin.xlsx",sheetName="intrinsic_subtye_population_only",append=T)
z.design <- matrix(c(
  c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
  c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
  c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
),ncol=5)
colnames(z.design) <- c("Luminial A","Luminal B",
                        "Luminal B HER2-",
                        "HER2 Enriched",
                        "Triple Negative")


parity <- data$parity_cat
x.covar <- cbind(parity,ethnicity.mat,refage)
model.lin.intrin <- EMmvpolySelfDesign(y.pheno.mis,x.self.design = x.covar[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.covar[,2:ncol(x.covar)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
model.lin.intrin[[4]][,2] <- c("Luminial A","Luminal B",
                               "Luminal B HER2-",
                               "HER2 Enriched",
                               "Triple Negative")
write.xlsx(model.lin.intrin[[4]],file="./risk_factor/result/risk_factor_paraty_lin.xlsx",sheetName="intrinsic_subtype_parity_lin",append=T)
#write.xlsx(model.lin[[5]],file="./risk_factor/result/risk_factor_paraty_lin.xlsx",sheetName="test_result",append=T)




z.design.new <- matrix(c(1-c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
                         c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)),ncol=2)
model.tri <- EMmvpolySelfDesign(y.pheno.mis,x.self.design = x.covar[,1,drop=F],z.design=z.design.new,baselineonly = NULL,additive = x.covar[,2:ncol(x.covar)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
model[[4]][,2] <- c("non-triple","triple")
write.xlsx(model.tri[[4]],file="./risk_factor/result/risk_factor_paraty_lin.xlsx",sheetName="non_triple_vs_triple",append=T)



# 
# y = y.pheno.mis
# additive = x.covar


# data.com <- fread("./data/concept_542-543-change-claude_bcac_pheno_v10_SPT_100217.txt")
# 
# 
# data.com.new <- data.com[,c(2,40,48,56,10)]
# 
# 
# new.
# 
# 
# getZstandard()