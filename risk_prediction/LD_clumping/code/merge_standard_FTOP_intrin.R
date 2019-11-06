setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/')
load('./intrinsic_subtypes_whole_genome/ICOG/result/meta_result_shared_1p.Rdata')
meta_intrin <- meta_result_shared_1p
load('./FTOP_whole_genome/ICOG/result/meta_result_shared_1p.Rdata')
meta_FTOP <- meta_result_shared_1p_FTOP
load("./standard_whole_genome/ICOG/result/meta_result_shared_1p.Rdata")

load("./EB_whole_genome/result/meta_result_shared_1p.Rdata")
#meta_stan <- meta_result_shared_1p
all.equal(meta_intrin$rs_id,meta_FTOP$rs_id)
all.equal(meta_stan$rs_id,meta_FTOP$rs_id)
meta_infor <- meta_intrin[,c(1:10,41:44)]
stan_result <- meta_stan[,c(11,12,17)]
colnames(stan_result) <- paste0("stan_",
                                c("logodds",
                                  "var",
                                  "p"))
FTOP_result <- meta_FTOP[,15]
colnames(meta_FTOP) <- c("FTOP")
intrin_result <- meta_intrin[,c(11:40,45)]
colnames(intrin_result)[1:5] <- c("Luminal_A",
                                  "Luminal_B",
                                  "Luminal_B_HER2Neg",
                                  "HER2_Enriched",
                                  "TN")
colnames(intrin_result)[31] <- "in_p"
# eb_result <- meta_eb[,c(11:16)]
# colnames(eb_result)[1:5] <- paste0("eb_",
#                                    c("Luminial_A",
#                                      "Luminal_B",
#                                      "Luminal_B_HER2Neg",
#                                      "HER2_Enriched",
#                                      "TN"))
# colnames(eb_result)[6] <- "heter_var"
whole_genome <- cbind(meta_infor,
                      stan_result,
                      FTOP_result,
                      intrin_result)
all.temp <- strsplit(whole_genome$var_name,split="_")
reference_allele <- rep("c",nrow(whole_genome))
effect_allele <- rep("c",nrow(whole_genome))

for(i in 1:nrow(whole_genome)){
  reference_allele[i] <- all.temp[[i]][[3]]
  effect_allele[i] <- all.temp[[i]][[4]]
}
whole_genome <- cbind(whole_genome,
                      reference_allele,
                      effect_allele)
whole_genome = whole_genome %>% 
  mutate(p.min = pmin(stan_p,FTOP_result,na.rm = T))
colnames(whole_genome)[19] <- "Luminal_A"
save(whole_genome,file = "./intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")
