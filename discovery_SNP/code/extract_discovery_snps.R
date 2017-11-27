library(R.utils)
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
x.covar.mis1 <- data1[,c(5:14,204)]
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
SG_ID <- SG_ID[idx.complete]




idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)
#load("./whole_genome_age/ICOG/ERPRHER2_fixed/result/score.test.support.icog.ERPRHER2.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")


snp.icogs.extract.id <- extract.list[,(ncol(extract.list)-2),drop=F]
write.table(snp.icogs.extract.id,file = paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_id_icog.txt"),quote = F,row.names=F)

Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed/"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)



qctool.command <- rep("c",564)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)


for(i in 1:564){
  geno.file <- Files[i]
  temp <- paste0("/spin1/users/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_id_icog.txt -og /spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2_extract",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/code/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)























#########fixed effect model
library(bc2)
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/LD_pruning.result")
extract.result.fixed <- extract.list.ld.pruning.result
idx <- which(is.na(extract.result.fixed$SNP.ICOGS)==T)

extract.result.fixed.shared <- extract.result.fixed[-idx,]
extract.result.onco.only <- extract.result.fixed[idx,]


load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_result_shared.Rdata")
extract.result.onco.fixed <- extract.result

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_result_shared.Rdata")
extract.result.icog.fixed <- extract.result



fill.match <- FillandMatch(extract.result.fixed.shared$SNP.ICOGS,extract.result.icog.fixed[[1]])
idx.fil <- fill.match[[1]]
idx.match <- fill.match[[2]]
icog.name <- extract.result.icog.fixed[[1]][idx.fil][idx.match]
all.equal(icog.name,extract.result.fixed.shared$SNP.ICOGS)
discover.data.icog.fixed <- extract.result.icog.fixed[[2]][,idx.fil][,idx.match]
colnames(discover.data.icog.fixed) <- icog.name
write.csv(discover.data.icog.fixed,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/discover.data.icog.fixed.csv")

save(discover.data.icog.fixed,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/discover.data.icog.fixed.Rdata")


fill.match <- FillandMatch(extract.result.fixed.shared$SNP.ONCO,extract.result.onco.fixed[[1]])
idx.fil <- fill.match[[1]]
idx.match <- fill.match[[2]]
onco.name <- extract.result.onco.fixed[[1]][idx.fil][idx.match]
all.equal(onco.name,extract.result.fixed.shared$SNP.ONCO)
discover.data.onco.fixed <- extract.result.onco.fixed[[2]][,idx.fil][,idx.match]
colnames(discover.data.onco.fixed) <- onco.name
write.csv(discover.data.onco.fixed,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/discover.data.onco.fixed.csv")

save(discover.data.onco.fixed,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/discover.data.onco.fixed.Rdata")






fill.match <- FillandMatch(extract.result.onco.only$SNP.ONCO,extract.result.onco.fixed[[1]])
idx.fil <- fill.match[[1]]
idx.match <- fill.match[[2]]
onco.name <- extract.result.onco.fixed[[1]][idx.fil][idx.match]
all.equal(onco.name,extract.result.onco.only$SNP.ONCO)
discover.data.onco.fixed.oncoonly <- extract.result.onco.fixed[[2]][idx.fil][idx.match]
colnames(discover.data.onco.fixed.oncoonly) <- onco.name
write.csv(discover.data.onco.fixed.oncoonly,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/discover.data.onco.fixed.oncoonly.csv")

save(discover.data.onco.fixed.oncoonly,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/discover.data.onco.fixed.oncoonly.Rdata")







idx <- which(extract.result.icog.fixed[[1]]%in%extract.result.fixed$SNP.ICOGS)
icog.snp.name <- as.character(extract.result.icog.fixed[[1]])[idx]
idx.match <- match(extract.result.fixed.shared$SNP.ICOGS,icog.snp.name)
icog.snp.name <- 
  icog.data.fixed <- extract.result.icog.fixed[[2]][,idx]

all.equal(icog.snp.name,extract.result.fixed.shared$SNP.ICOGS)
cbind(icog.snp.name,extract.result.fixed.shared$SNP.ICOGS)























load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/LD_pruning_casecase.result")
extract.result.random <- extract.list.ld.pruning.result


idx <- which(extract.result.fixed$p.value<=5E-08)
length(idx)
jdx <- which(extract.result.random$p.value<=5E-08)
length(jdx)

var.name.fixed <- extract.result.fixed$var_name
var.name.random <- extract.result.random$var_name
var.name.comb <- c(as.character(var.name.fixed),as.character(var.name.random))
var.name.comb <- unique(var.name.comb)
length(var.name.fixed)
length(var.name.random)
length(var.name.comb)


