#Goal: Generate all mahatattan and qq plots for the discovery paper
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
print(i1)
print(i2)
library(qqman)
library(data.table)
library(dplyr)

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
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
if(i1 ==1){
  
  data <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))
  gwas_result <- data %>%
    # filter(p.meta<=10^-4) %>% 
    select(c(SNP.Onco,chr.Onco,Position.Onco,p.meta))
  colnames(gwas_result) <- c("SNP","CHR",
                             "BP","P")
  gwas_result = gwas_result %>%
    mutate(P= ifelse(P==0,1E-20,P))
  
  fine_mapping <- read.csv("./data/filter_regions_standard.csv",header= T)
  gwas_result_filter <- FilterSNP(gwas_result,fine_mapping) 
  # dim(gwas_result)
  #idx <- sample(c(1:10760767),100000)
  if(i2==1){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_standard.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(gwas_result, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F)
    dev.off()
    
  }else if(i2==2){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_standard_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(gwas_result_filter, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F,ylim  = c(0,12))
    dev.off()
  }else if(i2==3){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_standard.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(gwas_result$P)
    dev.off()
  }else if(i2==4){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_standard_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(gwas_result_filter$P,ylim  = c(0,12))
    dev.off()
  }
  
}else if(i1==2){
  #load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
  meta_result_shared_1p <- as.data.frame(fread("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_fixed.txt",header=T))
  ftop <- meta_result_shared_1p
  ftop.p <- ftop$p.value
  #load("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p.Rdata")
  meta_result_shared_1p <- as.data.frame(fread("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_mixed.txt",header=T))
  mtop <- meta_result_shared_1p
  mtop.p <- mtop$p.value
  subtypes.p <- apply(cbind(ftop.p,mtop.p),1,min)
  idx <- which(is.na(ftop.p))
  ftop$subtypes.p <- subtypes.p
  ftop =  ftop %>% 
    mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0)|is.na(subtypes.p),1E-20,subtypes.p))
  subtypes_gwas_result <- ftop %>%
    select(rs_id,CHR,position,subtypes.p.new)
  colnames(subtypes_gwas_result) <- c("SNP",
                                      "CHR",
                                      "BP",
                                      "P")
  fine_mapping <- read.csv("./data/filter_regions_subtypes.csv",header= T)
  subtypes_gwas_result_filter <- FilterSNP(subtypes_gwas_result,fine_mapping) 
  
  if(i2==1){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_subtypes.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(subtypes_gwas_result, col= c("blue4", "orange3",),
              cex = 0.6,suggestiveline = F)
    dev.off()
    
  }else if(i2==2){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_subtypes_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(subtypes_gwas_result_filter, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F,ylim  = c(0,12))
    dev.off()
  }else if(i2==3){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_subtypes.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(subtypes_gwas_result$P)
    dev.off()
  }else if(i2==4){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_subtypes_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(subtypes_gwas_result_filter$P,ylim  = c(0,12))
    dev.off()
  }
}else if(i1==3){
  cimba_result_all <- as.data.frame(fread("./data/CIMBA_BCAC_meta_analysis_083019.txt",header = T))
  # idx <- which(cimba_result_all$CHR==11&
  #                cimba_result_all$position==132959475)
  # cimba_result_all[7617599:7617603,]
  # temp <- strsplit(cimba_result_all$MarkerName,
  #                  "_")
  # #temp[[1]][1]
  # n.snp <- nrow(cimba_result_all)
  # CHR <- rep(0,n.snp)
  # pos <- rep(0,n.snp)
  # for(i in 1:n.snp){
  #   if(i%%10000==0){
  #     print(i)
  #   }
  #   CHR[i] <- as.numeric(temp[[i]][1])
  #   pos[i] <- as.numeric(temp[[i]][2])
  # }
  # cimba_result_all$CHR <- CHR
  # cimba_result_all$position <- pos
  # cimba_result_all <- as.data.frame(cimba_result)
  
  
  colnames(cimba_result_all)[10] <- "P"
  cimba_result = cimba_result_all %>% 
    filter(Freq1>=0.008&
             Freq1<=0.992) %>% 
    select(MarkerName,CHR,
           position,P)
  
  colnames(cimba_result) <- c("SNP",
                              "CHR",
                              "BP",
                              "P")
  fine_mapping <- read.csv("./data/filter_regions_cimba.csv",header= T)
  cimba_result_filter <- FilterSNP(cimba_result,fine_mapping) 
  if(i2==1){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_cimba.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(cimba_result, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F,ylim  = c(0,300))
    dev.off()
    
  }else if(i2==2){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_cimba_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(cimba_result_filter, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F,ylim  = c(0,12))
    dev.off()
  }else if(i2==3){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_cimba.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(cimba_result$P,ylim  = c(0,300))
    dev.off()
  }else if(i2==4){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_cimba_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(cimba_result_filter$P,ylim  = c(0,12))
    dev.off()
  }
}else if(i1==4){
  #load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
  meta_result_shared_1p <- as.data.frame(fread("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_fixed.txt",header=T))
  ftop <- meta_result_shared_1p
  ftop.p <- ftop$p.value
  subtypes.p <- ftop.p
  ftop$subtypes.p <- ftop$p.value
  ftop =  ftop %>% 
    mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0)|is.na(subtypes.p),1E-20,subtypes.p))
  subtypes_gwas_result <- ftop %>%
    select(rs_id,CHR,position,subtypes.p.new)
  colnames(subtypes_gwas_result) <- c("SNP",
                                      "CHR",
                                      "BP",
                                      "P")
  fine_mapping <- read.csv("./data/filter_regions_subtypes.csv",header= T)
  subtypes_gwas_result_filter <- FilterSNP(subtypes_gwas_result,fine_mapping) 
  
  if(i2==1){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_fixed_subtypes.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(subtypes_gwas_result, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F)
    dev.off()
    
  }else if(i2==2){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_fixed_subtypes_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(subtypes_gwas_result_filter, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F,ylim  = c(0,12))
    dev.off()
  }else if(i2==3){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_fixed_subtypes.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(subtypes_gwas_result$P)
    dev.off()
  }else if(i2==4){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_fixed_subtypes_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(subtypes_gwas_result_filter$P,ylim  = c(0,12))
    dev.off()
  }
}else if(i1==4){
  #load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
  meta_result_shared_1p <- as.data.frame(fread("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_mixed.txt",header=T))
  ftop <- meta_result_shared_1p
  ftop.p <- ftop$p.value
  subtypes.p <- ftop.p
  ftop$subtypes.p <- ftop$p.value
  ftop =  ftop %>% 
    mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0),1E-20,subtypes.p))
  subtypes_gwas_result <- ftop %>%
    select(rs_id,CHR,position,subtypes.p.new)
  colnames(subtypes_gwas_result) <- c("SNP",
                                      "CHR",
                                      "BP",
                                      "P")
  fine_mapping <- read.csv("./data/filter_regions_subtypes.csv",header= T)
  subtypes_gwas_result_filter <- FilterSNP(subtypes_gwas_result,fine_mapping) 
  
  if(i2==1){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_fixed_subtypes.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(subtypes_gwas_result, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F)
    dev.off()
    
  }else if(i2==2){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_fixed_subtypes_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(subtypes_gwas_result_filter, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F,ylim  = c(0,12))
    dev.off()
  }else if(i2==3){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_fixed_subtypes.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(subtypes_gwas_result$P)
    dev.off()
  }else if(i2==4){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_fixed_subtypes_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(subtypes_gwas_result_filter$P,ylim  = c(0,12))
    dev.off()
  }
}else if(i1==5){
  #load("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p.Rdata")
  meta_result_shared_1p <- as.data.frame(fread("./whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_mixed.txt",header=T))
  mtop <- meta_result_shared_1p
  mtop.p <- mtop$p.value
  #just for easy coding, put mtop as ftop
  ftop = mtop
  ftop$subtypes.p <- ftop$p.value
  ftop =  ftop %>% 
    mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0),1E-20,subtypes.p))
  subtypes_gwas_result <- ftop %>%
    select(rs_id,CHR,position,subtypes.p.new)
  colnames(subtypes_gwas_result) <- c("SNP",
                                      "CHR",
                                      "BP",
                                      "P")
  fine_mapping <- read.csv("./data/filter_regions_subtypes.csv",header= T)
  subtypes_gwas_result_filter <- FilterSNP(subtypes_gwas_result,fine_mapping) 
  
  if(i2==1){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_mixed_subtypes.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(subtypes_gwas_result, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F)
    dev.off()
    
  }else if(i2==2){
    png(paste0("./discovery_SNP/result/manhattan_plot/man_mixed_subtypes_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    manhattan(subtypes_gwas_result_filter, col= c("blue4", "orange3"),
              cex = 0.6,suggestiveline = F,ylim  = c(0,12))
    dev.off()
  }else if(i2==3){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_mixed_subtypes.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(subtypes_gwas_result$P)
    dev.off()
  }else if(i2==4){
    png(paste0("./discovery_SNP/result/manhattan_plot/qq_mixed_subtypes_filter.png"),width = 7.635,height =4.7175,units = "in",res = 300,type="cairo")
    qq(subtypes_gwas_result_filter$P,ylim  = c(0,12))
    dev.off()
  }
}

