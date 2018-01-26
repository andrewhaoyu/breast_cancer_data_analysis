library(data.table)
setwd('/spin1/users/zhangh24/breast_cancer_data_analysis/')
hap3 <- fread("./data/HM3snplist_withCHRBP.txt",header=T)
dim(hap3)
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")



hap3.chr.bp <- cbind(hap3[,4],hap3[,5])

hap3.chr.bp.uni <- apply(hap3.chr.bp,1,paste0,collapse = ":")
dim(hap3)
length(hap3.chr.bp.uni)
hap3$chr.pos <- hap3.chr.bp.uni


bcac.chr.bp <- cbind(meta_result_shared_1p$CHR,
                     meta_result_shared_1p$position)
bcac.chr.bp.uni <- apply(bcac.chr.bp,1,paste0,collapse=":")

bcac <- cbind(meta_result_shared_1p$SNP.ICOGS,
              meta_result_shared_1p$SNP.ONCO,
              meta_result_shared_1p$CHR,
              meta_result_shared_1p$position,
              bcac.chr.bp.uni)












# share.id <- intersect(bcac.chr.bp.uni,hap3.chr.bp.uni)
# 
# hap3.in.bcac <- hap3.chr.bp.uni%in%share.id
# hap3.shared.uni <- hap3.chr.bp.uni[hap3.in.bcac]
# 
# bcac.in.hap3 <- bcac.chr.bp.uni%in%share.id
# bcac.shared.uni <- bcac.chr.bp.uni[bcac.in.hap3]
# 
# 
# 

extract.list.icog <- meta_result_shared_1p$SNP.ICOGS[hap3.in.bcac]
