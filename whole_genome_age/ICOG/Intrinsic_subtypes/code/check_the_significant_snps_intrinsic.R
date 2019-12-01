#check the number of significant snps in intrinsic subtypes analysis
setwd("/data/zhangh24/breast_cancer_data_analysis/")
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_final.Rdata")
intrinsic_result <- meta_result_shared_1p
#take out all the SNPs within the fine-mapping regions
FilterSNP <- function(gwas_result,fine_mapping){
  idx_cut <- NULL
  start <- fine_mapping$start
  end <- fine_mapping$end
  CHR <- fine_mapping$CHR
  
  #fine_mapping
  for(i in 1:nrow(fine_mapping)){
    print(i)
    chr_temp <- CHR[i]
    start_temp <- start[i]
    end_temp <- end[i]
    idx <- which(gwas_result$CHR==chr_temp&gwas_result$BP>=start_temp&gwas_result$BP<=end_temp)
    idx_cut <- c(idx_cut,idx)
  }
  ############duplicate variables remove by unqiue function
  idx_cut <- unique(idx_cut)
  gwas_result_filter <- gwas_result[-idx_cut,]
  return(gwas_result_filter)
}
library(dplyr)
intrinsic_gwas_result <- intrinsic_result %>%
  select(rs_id,CHR,position,p.value)
colnames(intrinsic_gwas_result) <- c("SNP",
                                    "CHR",
                                    "BP",
                                    "P")
fine_mapping <- read.csv("./data/filter_region_intrinsic_subtypes_analysis.csv",header= T)
colnames(fine_mapping)[3] <- "CHR"
intrinsic_gwas_result_sig <- intrinsic_gwas_result %>% filter(intrinsic_gwas_result$P<=5E-08)
intrinsic_gwas_result_filter <- FilterSNP(intrinsic_gwas_result_sig,fine_mapping) 
gwas_result <- intrinsic_gwas_result_filter
#Pick and keep the genome-wide significant SNPs within 500kb
idx.order <- order(gwas_result$CHR,gwas_result$BP)
gwas_result <- gwas_result[idx.order,]
PickUp <- function(gwas_result){
  result <- NULL
  #temporary top results within a region
  #if another top signal comes, then update with the new one
  #greedy algorithm
  temp.top <- NULL
  temp.top <- gwas_result[1,]
  for(i in 1:nrow(gwas_result)){
    temp <- gwas_result[i,]
    #if it is a SNP with lower P-value within the +-500kb region of the top SNP, then this SNP replace the top SNP
    if(temp$P<temp.top$P&
       temp$CHR==temp.top$CHR&
       temp$BP>=(temp.top$BP-500000)&
       temp$BP<=(temp.top$BP+500000)){
      #print(i)
     # break
      temp.top = temp
    }else if(temp$CHR!=temp.top$CHR|
             (temp$BP<=(temp.top$BP-500000))
             |
             (temp$BP>=(temp.top$BP+500000))){
      #if the SNP is out of the nearby region of the top snp, then the algorithm restarts
     # print(i)
     # break
      result <- rbind(result,temp.top)
      temp.top <- temp
    }
  }
return(result)  
}
significant_snp <- PickUp(gwas_result)

#load the subtyps analysis results
meta_result_shared_1p <- as.data.frame(fread("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_fixed.txt",header=T))
ftop <- meta_result_shared_1p
ftop.p <- ftop$p.value

meta_result_shared_1p <- as.data.frame(fread("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_mixed.txt",header=T))
mtop <- meta_result_shared_1p
mtop.p <- mtop$p.value
subtypes.p <- apply(cbind(ftop.p,mtop.p),1,min)
ftop$subtypes.p <- subtypes.p
ftop$mtop.p <- mtop$p.value
colnames(ftop)[1] <- "SNP"
snp_all = left_join(significant_snp,ftop,by="SNP")
head(snp_all)
snp_all = snp_all %>% select(SNP,CHR.x,BP,P,p.value,mtop.p,subtypes.p) %>% rename(CHR=CHR.x,intrinsic_P = P,ftop.p=p.value)

snp <- read.csv("./data/210_known_discovery_snp_paper_order.csv")
#take the subtypes snp
snp <- snp[201:208,]

snp = snp %>% 
  mutate(chr.pos = paste0(CHR,":",position))
intrinsic_gwas_result = intrinsic_gwas_result %>% 
  mutate(chr.pos = paste0(CHR,":",BP))
snp_dis = left_join(snp,intrinsic_gwas_result,
                    by="chr.pos") %>% 
  select(SNP,P)
head(snp_dis)
snp_dis_all <- left_join(snp_dis,ftop,by="SNP")
snp_dis_all = snp_dis_all %>% 
  select(SNP,CHR,position,P,p.value,mtop.p,subtypes.p) %>% 
  rename(BP=position,intrinsic_P=P,ftop.p=p.value)
snp_result <- rbind(snp_all,snp_dis_all)
write.csv(snp_result,file = paste0("./whole_genome_age/ICOG/Intrinsic_subtypes/result/significant_snp_intrinsic_subtypes.csv"))
