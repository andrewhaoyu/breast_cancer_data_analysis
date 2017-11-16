#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2_fixed/result/extract_result.Rdata")

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")


library(data.table)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
x.test.all.mis2 <- data2[,c(27:203,205)]

idx.control <- which(data2$Behaviour1==0)

x.test.all.mis2.control <- x.test.all.mis2[idx.control,]

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_result_shared.Rdata")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")

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

for(i in 1:ncol(extract.result.onco.control)){
  print(i)
  temp <- extract.result.onco.control[,i]
  
  for(j in 1:178){
    check.ld <- cor(x.test.all.mis2.control[,j],temp)^2
    if(check.ld>=0.1){
      idx.known.ld.flag <- c(idx.known.ld.flag,i)
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






# idx.extract.id <- c(idx.known.ld.flag,idx.julie.ld.flag)
# idx.extract.id <- unique(idx.extract.id)


extract.list.ld <- extract.list[idx.known.ld.flag,]

save(extract.list.ld,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract.list.ld")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract.list.ld")

save(idx.known.ld.flag,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/idx.known.ld.flag")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/idx.known.ld.flag")



LD.matrix <- cor(extract.result.onco.control)^2

#save(LD.matrix,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/LD.matrx.Rdata")

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/LD.matrx.Rdata")



LD_pruning = function(sig_SNPs,LD2){
  sig_SNPs_temp =sig_SNPs
  filter_result = NULL
  LD2.temp =LD2
  temp.ind = 1
  while(nrow(sig_SNPs_temp)!=0){
    
    idx = which.min(sig_SNPs_temp$p.value)
    filter_result = rbind(filter_result,sig_SNPs_temp[idx,])
    LD2.single = LD2.temp[idx,]
    idx.cut = which( LD2.single>=0.1)
    position.range <- 500*10^3
    filter_result_position = sig_SNPs_temp$position[idx]
    filter_CHR = sig_SNPs_temp$CHR[idx]
    idx.cut2 <- which((sig_SNPs_temp$position>=filter_result_position-position.range)&(sig_SNPs_temp$position<=filter_result_position+position.range)&(sig_SNPs_temp$CHR==filter_CHR ))
    idx.cut <- c(idx.cut,idx.cut2)
    idx.cut <- unique(idx.cut)
    LD2.temp = LD2.temp[-idx.cut,-idx.cut]
    LD2.temp = as.matrix(LD2.temp)
    sig_SNPs_temp = sig_SNPs_temp[-idx.cut,]
    temp.ind = temp.ind+1
  }
  return(filter_result)
}



LD.matrix <- LD.matrix[-idx.known.ld.flag,-idx.known.ld.flag]
extract.list <- extract.list[-idx.known.ld.flag,]















extract.list.ld.pruning.result <- LD_pruning(extract.list,LD.matrix)

save(extract.list.ld.pruning.result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/LD_pruning.result")

write.csv(extract.list.ld.pruning.result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/LD_pruning.csv",quote=F)










