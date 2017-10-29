#load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_result.Rdata")

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/extract_result.Rdata")

data2 <- read.csv("./data/ONCO_pruning.csv",header=T)
x.test.all.mis2 <- data2[,c(27:205)]

idx.control <- which(data2$Behaviour1==0)

x.test.all.mis2.control <- x.test.all.mis2[idx.control,]

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
# extract.list <- extract.list[-c(1700,1701,1702,1703,1705,1706,1707),]
# save(extract.list,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
#try <- duplicated(extract.list)

extract.result.onco <- extract.result[[2]]
idx.match <- match(extract.list$SNP.ONCO,extract.result[[1]])
extract.result.onco <- extract.result.onco[,idx.match]
extract.result.onco.control <- extract.result.onco[idx.control,]


idx.known.ld.flag <- NULL

for(i in 1:ncol(extract.result.onco.control)){
  print(i)
  temp <- extract.result.onco.control[,i]
  
  for(j in 1:179){
    check.ld <- cor(x.test.all.mis2.control[,j],temp)^2
    if(check.ld>=0.1){
      idx.known.ld.flag <- c(idx.known.ld.flag,i)
      break
    }
      
  }
  
  
  
}

extract.list.ld <- extract.list[idx.known.ld.flag,]

save(extract.list.ld,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract.list.ld")








LD.matrix <- cor(extract.result.onco.control)^2

save(LD.matrix,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/LD.matrx.Rdata")



LD_pruning = function(sig_SNPs,LD2){
  sig_SNPs_temp =sig_SNPs
  filter_result = NULL
  LD2.temp =LD2
  temp.ind = 1
  while(nrow(sig_SNPs_temp)!=0){
    
    idx = which.min(sig_SNPs_temp$pvalue)
    filter_result = rbind(filter_result,sig_SNPs_temp[idx,])
    LD2.single = LD2.temp[idx,]
    idx.cut = which( LD2.single>=0.1)
    position.range <- 10^6
    filter_result_position = sig_SNPs_temp$position[idx]
    idx.cut2 <- which((sig_SNPs_temp$position>=filter_result_position-position.range)&(sig_SNPs_temp$position<=filter_result_position+position.range))
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



new_filter <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))

idx_cut <- NULL
position.cut <- 10^6
for(i in 1:nrow(new_filter)){
  print(i)
  chr_temp <- new_filter[i,3]
  position_temp <- new_filter[i,2]
  position_low <- position_temp-position.cut
  position_high <- position_temp+position.cut
  idx <- which(extract.list$CHR==chr_temp&extract.list$position>position_low&
                 extract.list$position<position_high)
  idx_cut <- c(idx_cut,idx)
}
idx_cut <- unique(idx_cut)
extract.list <- extract.list[-idx_cut,]

LD.matrix <- LD.matrix[-idx_cut,-idx_cut]




extract.list <- LD_pruning(extract.list,LD.matrix)


save(extract.list,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/LD_pruning.result")



load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/LD_pruning.result")
top <- extract.list[extract.list$pvalue<=1e-07,]


load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_result.Rdata")
extract.result.icog <- extract.result[[2]]
extract.result.icog.id <- extract.result[[1]]

idx <- which((extract.result.icog.id%in%top$SNP.ICOGS)==T)
extract.result.icog.top <- extract.result.icog[,idx]
extract.result.icog.id.top <- extract.result.icog.id[idx]

load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/extract_result.Rdata")
extract.result.onco <- extract.result[[2]]
extract.result.onco.id <- extract.result[[1]]

idx <- which((extract.result.onco.id%in%top$SNP.ONCO)==T)
extract.result.onco.top <- extract.result.onco[,idx]
extract.result.onco.id.top <- extract.result.onco.id[idx]

cbind(extract.result.onco.id.top,extract.result.onco.id.top)

colnames(extract.result.icog.top) <- extract.result.icog.id.top
colnames(extract.result.onco.top) <- 
  extract.result.onco.id.top
save(extract.result.icog.top,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/topsnpvalue.Rdata")
save(extract.result.onco.top,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/topsnpvalue.Rdata")








load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_result.Rdata")

