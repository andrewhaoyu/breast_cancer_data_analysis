#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2_fixed/result/extract_result.Rdata")

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")


library(data.table)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
x.test.all.mis2 <- data2[,c(27:203,205)]

idx.control <- which(data2$Behaviour1==0)

x.test.all.mis2.control <- x.test.all.mis2[idx.control,]

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_casecase/result/extract_result_shared.Rdata")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_list.Rdata")

# extract.list <- extract.list[-c(1700,1701,1702,1703,1705,1706,1707),]
# save(extract.list,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
#try <- duplicated(extract.list)

extract.result.onco <- extract.result[[2]]
extract.result.onco.control <- extract.result.onco[idx.control,]

library(data.table)
Julie_snp <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_onco.csv"))
x.test.Julie <- Julie_snp[,-1]
x.test.Julie.control <- x.test.Julie[idx.control,]

idx.julie.ld.flag <- NULL


idx.known.ld.flag <- NULL
ld.known.id <- NULL

for(i in 1:ncol(extract.result.onco.control)){
  print(i)
  temp <- extract.result.onco.control[,i]
  
  for(j in 1:178){
    check.ld <- cor(x.test.all.mis2.control[,j],temp)^2
    if(check.ld>=0.1){
      idx.known.ld.flag <- c(idx.known.ld.flag,i)
      ld.known.id <- c(ld.known.id,j) 
      break
    }
    
  }
  for(j in 1:ncol(x.test.Julie)){
    check.ld <- cor(x.test.Julie.control[,j],temp)^2
    if(check.ld>=0.1){
      idx.julie.ld.flag <- c(idx.julie.ld.flag,i)
      break
    }
  }
  
  
  
}

ld.known.data <- data.frame(idx.known.ld.flag,ld.known.id,extract.list[idx.known.ld.flag,13:14],stringsAsFactors = F)

head(ld.known.data)













load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata")

fine_mapping <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/fine_mapping_annotated_clean.csv",header= T,
                         stringsAsFactors = F)

idx_cut <- NULL
start <- fine_mapping$start
end <- fine_mapping$end
CHR <- fine_mapping$CHR
position <- fine_mapping$position


# fine.info <- data.frame(CHR,position,
#                         stringsAsFactors = F)
# 
# 
# 
# 
# idx.fine <- which((do.call(paste0,all.info)%in%do.call(paste0,fine.info))==T)
# fine_mapping_snp_names <- meta_result_shared_1p[idx.fine,]


all <- meta_result_shared_1p





library(bc2)
idx.temp.known <- get_fine_mapping_id(meta_result_shared_1p,fine_mapping)
idx.cut <- idx.temp.known[[1]]

known.flag <- idx.temp.known[[2]]
all.known.region.snps <- meta_result_shared_1p[idx.cut,]
all.known.region.snps <- cbind(all.known.region.snps,known.flag)
library(dplyr)
dim(all.known.region.snps)
all.known.region.snps <- all.known.region.snps%>%filter(p.value<=0.01)


all.known.region.snps <- all.known.region.snps[,c(13,14,16)]

ld.known.data <- ld.known.data[,c(3,4,2)]
colnames(ld.known.data) <- colnames(all.known.region.snps)

all.known.region.snps <- rbind(all.known.region.snps,ld.known.data)


save(all.known.region.snps,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.known.region.snps.Rdata")

known.region.snps.icogs <- all.known.region.snps$SNP.ICOGS
known.region.snps.onco <- all.known.region.snps$SNP.ONCO
write.table(known.region.snps.icogs,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/known.region.icog.txt",quote=F,row.names=F)
write.table(known.region.snps.onco,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/known.region.onco.txt",quote=F,row.names=F)

discovery.info <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snps_annotated_clean.csv",header=T,stringsAsFactors = F)


idx.temp.dis <- get_fine_mapping_id(meta_result_shared_1p,discovery.info)
idx.cut <- idx.temp.dis[[1]]

dis.flag <- idx.temp.dis[[2]]
all.dis.region.snps <- meta_result_shared_1p[idx.cut,]
all.dis.region.snps <- cbind(all.dis.region.snps,dis.flag)
library(dplyr)
dim(all.dis.region.snps)
all.dis.region.snps <- all.dis.region.snps%>%filter(p.value <= 0.01)
dim(all.dis.region.snps)
all.dis.region.snps <- all.dis.region.snps[,c(13,14,16)]

save(all.known.region.snps,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.known.region.snps.Rdata")







dis.region.snps.icogs <- all.dis.region.snps$SNP.ICOGS
dis.region.snps.onco <- all.dis.region.snps$SNP.ONCO
write.table(dis.region.snps.icogs,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/dis.region.icog.txt",quote=F,row.names=F)
write.table(dis.region.snps.onco,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/dis.region.onco.txt",quote=F,row.names=F)




test <- which(all.dis.region.snps$SNP.ICOGS%in%all.known.region.snps$SNP.ICOGS==F)


test.1 <- which(all.dis.region.snps$SNP.ICOGS%in%all.known.region.snps$SNP.ICOGS==T)
length(test.1)


dis.region.snps <- all.dis.region.snps[test,]
dis.region.snps$dis.flag <- dis.region.snps$dis.flag+178

colnames(dis.region.snps) <- colnames(all.known.region.snps)
all.conditional.snps <- rbind(all.known.region.snps,dis.region.snps)
save(all.conditional.snps,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")

all.conditional.snps.icogs <- all.conditional.snps$SNP.ICOGS
all.conditional.snps.onco <- all.conditional.snps$SNP.ONCO
write.table(all.conditional.snps.icogs,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.icogs.txt",quote=F,row.names=F)
write.table(all.conditional.snps.onco,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.onco.txt",quote=F,row.names=F)





#which(is.na(new$SNP.ONCO)==T)

