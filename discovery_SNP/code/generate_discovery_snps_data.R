##############generate the discovery snps list

load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_result_shared.Rdata")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/LD_pruning.result")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")
discovery_snp <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary.csv",header=T)
discovery.icogs <- as.character(discovery_snp$SNP.ICOGS)
idx.fil <- which((extract.list[[1]]%in%discovery.icogs)==T)

idx.match <-match(discovery.icogs,extract.list[[1]][idx.fil])
























#########fixed effect model
library(bc2)
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/LD_pruning.result")
extract.result.fixed <- extract.list.ld.pruning.result
idx <- which(is.na(extract.result.fixed$SNP.ICOGS)==T)

extract.result.fixed.shared <- extract.result.fixed[-idx,]
extract.result.onco.only <- extract.result.fixed[idx,]


load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_result_shared.Rdata")
extract.result.onco.fixed <- extract.result

load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_result_shared.Rdata")
extract.result.icog.fixed <- extract.result



fill.match <- FillandMatch(extract.result.fixed.shared$SNP.ICOGS,extract.result.icog.fixed[[1]])
idx.fil <- fill.match[[1]]
idx.match <- fill.match[[2]]
icog.name <- extract.result.icog.fixed[[1]][idx.fil][idx.match]
all.equal(icog.name,extract.result.fixed.shared$SNP.ICOGS)
discover.data.icog.fixed <- extract.result.icog.fixed[[2]][,idx.fil][,idx.match]
colnames(discover.data.icog.fixed) <- icog.name
write.csv(discover.data.icog.fixed,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/discover.data.icog.fixed.csv")

save(discover.data.icog.fixed,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/discover.data.icog.fixed.Rdata")


fill.match <- FillandMatch(extract.result.fixed.shared$SNP.ONCO,extract.result.onco.fixed[[1]])
idx.fil <- fill.match[[1]]
idx.match <- fill.match[[2]]
onco.name <- extract.result.onco.fixed[[1]][idx.fil][idx.match]
all.equal(onco.name,extract.result.fixed.shared$SNP.ONCO)
discover.data.onco.fixed <- extract.result.onco.fixed[[2]][,idx.fil][,idx.match]
colnames(discover.data.onco.fixed) <- onco.name
write.csv(discover.data.onco.fixed,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/discover.data.onco.fixed.csv")

save(discover.data.onco.fixed,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/discover.data.onco.fixed.Rdata")






fill.match <- FillandMatch(extract.result.onco.only$SNP.ONCO,extract.result.onco.fixed[[1]])
idx.fil <- fill.match[[1]]
idx.match <- fill.match[[2]]
onco.name <- extract.result.onco.fixed[[1]][idx.fil][idx.match]
all.equal(onco.name,extract.result.onco.only$SNP.ONCO)
discover.data.onco.fixed.oncoonly <- extract.result.onco.fixed[[2]][idx.fil][idx.match]
colnames(discover.data.onco.fixed.oncoonly) <- onco.name
write.csv(discover.data.onco.fixed.oncoonly,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/discover.data.onco.fixed.oncoonly.csv")

save(discover.data.onco.fixed.oncoonly,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/discover.data.onco.fixed.oncoonly.Rdata")







idx <- which(extract.result.icog.fixed[[1]]%in%extract.result.fixed$SNP.ICOGS)
icog.snp.name <- as.character(extract.result.icog.fixed[[1]])[idx]
idx.match <- match(extract.result.fixed.shared$SNP.ICOGS,icog.snp.name)
icog.snp.name <- 
icog.data.fixed <- extract.result.icog.fixed[[2]][,idx]

all.equal(icog.snp.name,extract.result.fixed.shared$SNP.ICOGS)
cbind(icog.snp.name,extract.result.fixed.shared$SNP.ICOGS)























load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/LD_pruning_casecase.result")
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


