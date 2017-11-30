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

ld.known.data <- data.frame(idx.known.ld.flag,ld.known.id,extract.result[[1]][idx.known.ld.flag],stringsAsFactors = F)



load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata")

fine_mapping <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/fine_mapping_regions.csv",header= T)

idx_cut <- NULL
start <- fine_mapping$start
end <- fine_mapping$end
CHR <- fine_mapping$V3
position <- fine_mapping$V4

all.info <- data.frame(meta_result_shared_1p$CHR,
                  meta_result_shared_1p$position,
                  stringsAsFactors = F)
colnames(all.info) <- c("CHR","position")
fine.info <- data.frame(CHR,position,
                        stringsAsFactors = F)




idx.fine <- which((do.call(paste0,all.info)%in%do.call(paste0,fine.info))==T)
fine_mapping_snp_names <- meta_result_shared_1p[idx.fine,]

known.flag <- NULL


library(bc2)
idx.temp <- get_fine_mapping_id(meta_result_shared_1p,fine_mapping)

all.known.region.snps <- meta_result_shared_1p[idx_cut,]

library(dplyr)

 test <- all.known.region.snps%>%filter(p.value <= 1e-04)





new_filter <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/Filter_based_on_Montse.csv",header=T,stringsAsFactors = F)
new_filter[,2] <- as.numeric(gsub(",","",new_filter[,2]))

idx_cut <- NULL

position.cut <- 500*10^3

for(i in 1:nrow(new_filter)){
  print(i)
  chr_temp <- new_filter[i,3]
  position_temp <- new_filter[i,2]
  position_low <- position_temp-position.cut
  position_high <- position_temp+position.cut
  idx <- which(meta_result_shared_1p_filter$CHR==chr_temp&meta_result_shared_1p_filter$position>position_low&
                 meta_result_shared_1p_filter$position<position_high)
  idx_cut <- c(idx_cut,idx)
}
idx_cut <- unique(idx_cut)
meta_result_shared_1p_filter_Ju <- meta_result_shared_1p_filter[-idx_cut,]

save(meta_result_shared_1p_filter_Ju,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_filter_1M_Ju.Rdata")



































