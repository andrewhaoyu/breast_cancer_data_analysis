setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/risk_prediction/')
load('./intrinsic_subtypes_whole_genome/ICOG/result/meta_result_shared_1p.Rdata')
meta_intrin <- meta_result_shared_1p
load('./FTOP_whole_genome/ICOG/result/meta_result_shared_1p.Rdata')
meta_FTOP <- meta_result_shared_1p_FTOP
load("./standard_whole_genome/ICOG/result/meta_result_shared_1p.Rdata")
meta_stan <- meta_result_shared_1p
all.equal(meta_intrin$rs_id,meta_FTOP$rs_id)
all.equal(meta_stan$rs_id,meta_FTOP$rs_id)
meta_infor <- meta_intrin[,c(1:10,41:44)]
