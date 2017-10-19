setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
pheno.file <- "./data/pheno.onco"
load(pheno.file)
n.sub = nrow(pheno)
y.pheno.mis2 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1")


idx.fil <- onco.order[,1]%in%pheno$Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(pheno$Onc_ID,onco.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome/ONCO/ERPRHER2_fixed/result/score.test.support.onco.ERPRHER2.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
snp.onco.extract.id <- extract.list[,(ncol(extract.list)-1),drop=F]
write.table(snp.icogs.extract.id,file = paste0("./whole_genome/ONCO/ERPRHER2_fixed/result/extract_id_onco.txt"),quote = F,row.names=F)



Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
geno.file <- Files[i1]
