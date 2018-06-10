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
#idx.remove <- which(data$status==2|data$status==3)
#data <- data[-idx.remove,]
#########remove people with missing parity_cat


idx.remove <- which(data$parity_cat==9)
data <- data[-idx.remove,]


status <- data$status
status[(status==2|status==3)] = 1
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
write.xlsx(model.lin[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="second_stage_parameter")
write.xlsx(model.lin[[5]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="test_result",append=T)
# write.xlsx(model[[7]],file="./risk_factor/result/risk_factor_no_study.xlsx",sheetName="first_stage_result",append=T)



load("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/z.design.Rdata")

model.lin.tninter <- EMmvpolySelfDesign(y.pheno.mis,x.self.design = x.covar[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.covar[,2:ncol(x.covar)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

write.xlsx(model.lin.tninter[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="second_stage_parameter_third_order",append=T)
write.xlsx(model.lin.tninter[[5]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="test_result_third_order",append=T)

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
# design_cat <- data$design_cat
# #idx.pop <- which(data$design_cat==0)
# parity.mat_design_inter <- design_cat*parity.mat
# colnames(parity.mat_design_inter) <- paste0(colnames(parity.mat),"_interaction")
# x.covar <- cbind(parity.mat,parity.mat_design_inter,design_cat,ethnicity.mat,refage)
# model.intrinsic.design.inter = EMmvpolySelfDesign(y.pheno.mis,x.self.design = x.covar,z.design=z.design,
#                             baselineonly = NULL,additive =NULL,pairwise.interaction = NULL,
#                             saturated = NULL,
#                             missingTumorIndicator = 888)
# model2[[4]][,2] <- rep(c("Luminial A","Luminal B",
#                          "Luminal B HER2-",
#                          "HER2 Enriched",
#                          "Triple Negative"),15)
# write.xlsx(model2[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="intrinsic_subtye",append=T)

idx.pop <- which(data$design_cat==0)
# idx.non_missing <- which(data$molgroup!=9)
# data.pop.com <- data[idx.pop[idx.non_missing],]


# parity.mat_design_inter <- design_cat*parity.mat
# colnames(parity.mat_design_inter) <- paste0(colnames(parity.mat),"_interaction")
x.covar <- cbind(parity.mat,ethnicity.mat,refage)
model.intrinsic.pop = EMmvpolySelfDesign(y.pheno.mis[idx.pop,],x.self.design = x.covar[idx.pop,],z.design=z.design,
                            baselineonly = NULL,additive =NULL,pairwise.interaction = NULL,
                            saturated = NULL,
                            missingTumorIndicator = 888)
#model.intrinsic.pop <- model3
model.intrinsic.pop[[4]][,2] <- rep(c("Luminial A","Luminal B",
                         "Luminal B HER2-",
                         "HER2 Enriched",
                         "Triple Negative"),10)
write.xlsx(model.intrinsic.pop[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="intrinsic_subtye_population_only",append=T)

# z.design <- matrix(c(
#   c(0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0),
#   c(0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1),
#   c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
#   c(0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0),
#   c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)
# ),ncol=5)
# colnames(z.design) <- c("Luminial A","Luminal B",
#                         "Luminal B HER2-",
#                         "HER2 Enriched",
#                         "Triple Negative")
# 
# 
# parity <- data$parity_cat
# x.covar <- cbind(parity,ethnicity.mat,refage)
# model.lin.intrin <- EMmvpolySelfDesign(y.pheno.mis,x.self.design = x.covar[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.covar[,2:ncol(x.covar)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
# model.lin.intrin[[4]][,2] <- c("Luminial A","Luminal B",
#                                "Luminal B HER2-",
#                                "HER2 Enriched",
#                                "Triple Negative")
# write.xlsx(model.lin.intrin[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="intrinsic_subtype_parity_lin",append=T)
#write.xlsx(model.lin[[5]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="test_result",append=T)




# z.design.new <- matrix(c(1-c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
#                          c(1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)),ncol=2)
# model.tri <- EMmvpolySelfDesign(y.pheno.mis,x.self.design = x.covar[,1,drop=F],z.design=z.design.new,baselineonly = NULL,additive = x.covar[,2:ncol(x.covar)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
# model.tri[[4]][,2] <- c("non-triple","triple")
# write.xlsx(model.tri[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="non_triple_vs_triple",append=T)



###########################
data <- fread("./data/dataset_montse_20180524.txt")



idx.remove <- which(data$parity_cat==0|is.na(data$breastmos6)|
                      data$breastmos_cat==9)
data <- data[-idx.remove,]


status <- data$status
status[(status==2|status==3)] = 1
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
breast_mos.mat <- model.matrix(~as.factor(data$breastmos_cat)-1)[,-1]
colnames(breast_mos.mat) <- paste0("breast_cat",c(2:5))
study.mat <- model.matrix(~as.factor(data$study)-1)[,-1]
colnames(study.mat) <- unique(data$study)[-1]
ethnicity.mat <- model.matrix(~as.factor(data$ethnicity)-1)[,-1]
colnames(ethnicity.mat) <- paste0("ethnicity",c(1:6))[-1]
refage <- data$refage

breastmos6 <- data$breastmos6
breast_mos <- data$breastmos_cat
x.covar <- cbind( breast_mos,ethnicity.mat,refage)
model.lin.b <- TwoStageModel(y = y.pheno.mis,
                           #baselineonly = study.mat,
                           additive = x.covar,
                           missingTumorIndicator = 888)
write.xlsx(model.lin.b[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="second_stage_parameter_breast",append = T)
write.xlsx(model.lin.b[[5]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="test_result_breast",append=T)
#write.csv(data[idx,],file="./risk_factor/result/special_people.csv")

load("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/z.design.Rdata")

model.lin.tninter.b <- EMmvpolySelfDesign(y.pheno.mis,x.self.design = x.covar[,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.covar[,2:ncol(x.covar)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

write.xlsx(model.lin.tninter.b[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="sparameter_third_order_breast",append=T)
write.xlsx(model.lin.tninter.b[[5]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="test_result_third_order_breast",append=T)




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

idx.pop <- which(data$design_cat==0)

x.covar <- cbind(breast_mos.mat,ethnicity.mat,refage)
model.intrinsic.pop.b = EMmvpolySelfDesign(y.pheno.mis[idx.pop,],x.self.design = x.covar[idx.pop,],z.design=z.design,
                                         baselineonly = NULL,additive =NULL,pairwise.interaction = NULL,
                                         saturated = NULL,
                                         missingTumorIndicator = 888)
#model.intrinsic.pop <- model3
model.intrinsic.pop.b[[4]][,2] <- rep(c("Luminial A","Luminal B",
                                      "Luminal B HER2-",
                                      "HER2 Enriched",
                                      "Triple Negative"),10)
write.xlsx(model.intrinsic.pop.b[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="is_population_only_breast",append=T)
























data <- fread("./data/dataset_montse_20180524.txt")



idx.remove <- which(data$parity_cat==9)
data <- data[-idx.remove,]

idx.nob <- which(data$breastmos_cat==0)
idx.b <- which(data$breastmos_cat>=1&data$breastmos_cat<=5)
idx.m <- which(data$breastmos_cat==9)
idx.pb <- which(data$breastmos_cat>=1&data$breastmos_cat<=5&data$design_cat==0)
idx.pm <- which(data$breastmos_cat==9&data$design_cat==0)

status <- data$status
status[(status==2|status==3)] = 1
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
#breast_mos.mat <- model.matrix(~as.factor(data$breastmos_cat)-1)[,-(1:2)]
#colnames(breast_mos.mat) <- paste0("breast_cat",c(2:5))
parity_cat <- data$parity_cat
parity.mat <- model.matrix(~as.factor(data$parity_cat)-1)[,-1]
colnames(parity.mat) <- paste0("parity_cat",c(1:4))

study.mat <- model.matrix(~as.factor(data$study)-1)[,-1]
colnames(study.mat) <- unique(data$study)[-1]
ethnicity.mat <- model.matrix(~as.factor(data$ethnicity)-1)[,-1]
colnames(ethnicity.mat) <- paste0("ethnicity",c(1:6))[-1]
refage <- data$refage

library(xlsx)
#breast_mos <- data$breastmos_cat
x.covar <- cbind(parity_cat,ethnicity.mat,refage)

model.p.b <- TwoStageModel(y = y.pheno.mis[idx.b,],
                             #baselineonly = study.mat,
                             additive = x.covar[idx.b,],
                             missingTumorIndicator = 888)
write.xlsx(model.p.b[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="ssp_pwithb",append = T)
write.xlsx(model.p.b[[5]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="test_result_pwithb",append=T)

idx.rare <- which(ethnicity.mat[idx.m,1]==1|ethnicity.mat[idx.m,3]==1)
model.p.mb <- TwoStageModel(y = y.pheno.mis[idx.m,],
                           baselineonly = x.covar[idx.m,5,drop=F],
                           additive = x.covar[idx.m,
                                              c(1,2,3,4,6,7)],
                           missingTumorIndicator = 888)
write.xlsx(model.p.mb[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="ssp_pwithmb",append = T)
write.xlsx(model.p.mb[[5]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="test_result_pwithmb",append=T)

#write.csv(data[idx,],file="./risk_factor/result/special_people.csv")

load("./known_SNPs/known_SNPs_analysis_G_revised/additive_model_third_order/result/z.design.Rdata")

model.lin.tninter.p.b <- EMmvpolySelfDesign(y.pheno.mis[idx.b,],x.self.design = x.covar[idx.b,1,drop=F],z.design=z.design,baselineonly = NULL,additive = x.covar[idx.b,2:ncol(x.covar)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)

write.xlsx(model.lin.tninter.p.b[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="sparameter_third_order_pwithb",append=T)
write.xlsx(model.lin.tninter.p.b[[5]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="test_result_third_order_pwithb",append=T)

model.lin.tninter.p.mb <- EMmvpolySelfDesign(y.pheno.mis[idx.m,],x.self.design = x.covar[idx.m,1,drop=F],z.design=z.design,baselineonly = x.covar[idx.m,5,drop=F],additive = x.covar[idx.m,c(2,3,4,6,7)],pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)
write.xlsx(model.lin.tninter.p.mb[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="sparameter_third_order_pwithmb",append=T)
write.xlsx(model.lin.tninter.p.mb[[5]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="test_result_third_order_pwithmb",append=T)



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

idx.pop <- which(data$design_cat==0)
#idx.rare <- which(x.covar[idx.pb,5]==1|
                    #x.covar[idx.pb,7]==1)
x.covar <- cbind(parity.mat,ethnicity.mat,refage)
model.intrinsic.pop.pb = EMmvpolySelfDesign(y.pheno.mis[idx.pb,],x.self.design = x.covar[idx.pb,c(5:10)],z.design=z.design,
                                           baselineonly = NULL,additive = NULL,pairwise.interaction = NULL,
                                           saturated = NULL,
                                           missingTumorIndicator = 888)
#model.intrinsic.pop <- model3
model.intrinsic.pop.pb[[4]][,2] <- rep(c("Luminial A","Luminal B",
                                        "Luminal B HER2-",
                                        "HER2 Enriched",
                                        "Triple Negative"),7)
write.xlsx(model.intrinsic.pop.pb[[4]],file="./risk_factor/result/risk_factor_all_include.xlsx",sheetName="is_population_only_breast",append=T)

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