#prepare data for locus plot
setwd("/Users/haoyuzhang/GoogleDrive/breast_cancer_data_analysis")
#standard analysis data
data <- as.data.frame(fread("/Users/haoyuzhang/GoogleDrive/BCAC\ OncoArray\ Subtypes\ Local/breast_cancer_paper_writing/breast_cancer_discovery_paper/core_figures/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))
discovery_snp <- read.csv("./data/discovery_snp_summary_new.csv",header=T)
library(dplyr)
data_plot = data %>% select(
  c(var_name,SNP.Onco,chr.Onco,Position.Onco,p.meta
))
colnames(data_plot) = c("var_name",
                        "MarkerName",
                        "CHR",
                        "position",
                        "P.value")
temp.str <- strsplit(data_plot$MarkerName,":")
rs_id <- rep("c",nrow(data_plot))
for(i in 1:nrow(data_plot)){
  rs_id[i] <- temp.str[[i]][1]
}
data_plot$rs_id <- rs_id
data_plot <- data_plot %>% select(
  c(rs_id,CHR,position,P.value)
)
colnames(data_plot)[1] <- "MarkerName"
#write.table(data_plot,file = "./discovery_SNP/result/locus_plot_data/data_plot.txt",col.names=T,row.names=F,quote=F)
idx <- NULL
for(k in 1:nrow(discovery_snp)){
  idx.temp <- which(data_plot$CHR==discovery_snp$CHR.x[k]&
                 data_plot$position>= discovery_snp$Pos[k]-10^6&
                 data_plot$position<= discovery_snp$Pos[k]+10^6)
  data_plot_try <- data_plot[idx.temp,]
  #if the SNP don't have rs_id
  #use the chr:position format
  jdx <- which(c(1:nrow(data_plot_try)%in%grep("rs",data_plot_try$MarkerName))!=T)
  data_plot_try[jdx,1] <- paste0("chr",data_plot_try$CHR,":",data_plot_try$position)[jdx]
  write.table(data_plot_try,file = paste0("./discovery_SNP/result/locus_plot_data/data_plot_",discovery_snp$SNP.ONCO[k],".txt"),sep="\t",col.names=T,row.names=F,quote=F)
  
}
#take out SNPs in LD with known SNP
k <- 14
idx.temp <- which(data_plot$CHR==discovery_snp$CHR.x[k]&
                    data_plot$position>= discovery_snp$Pos[k]-10^6&
                    data_plot$position<= discovery_snp$Pos[k]+10^6)
data_plot_try <- data_plot[idx.temp,]
#if the SNP don't have rs_id
#use the chr:position format
jdx <- which(c(1:nrow(data_plot_try)%in%grep("rs",data_plot_try$MarkerName))!=T)
data_plot_try[jdx,1] <- paste0("chr",data_plot_try$CHR,":",data_plot_try$position)[jdx]
kdx <- which(data_plot_try$CHR==discovery_snp$CHR.x[k]&
               data_plot_try$position>= discovery_snp$Pos[k]-10^5&
               data_plot_try$position<= discovery_snp$Pos[k]+10^5&
               data_plot_try$P.value<1.49612e-10)
data_plot_try <- data_plot_try[-kdx,]
write.table(data_plot_try,file = paste0("./discovery_SNP/result/locus_plot_data/data_plot_",discovery_snp$SNP.ONCO[k],".txt"),sep="\t",col.names=T,row.names=F,quote=F)
# k <- 22
# idx.temp <- which(data_plot$CHR==discovery_snp$CHR.x[k]&
#                     data_plot$position>= discovery_snp$Pos[k]-10^6&
#                     data_plot$position<= discovery_snp$Pos[k]+10^6)
# data_plot_try <- data_plot[idx.temp,]
# #if the SNP don't have rs_id
# #use the chr:position format
# jdx <- which(c(1:nrow(data_plot_try)%in%grep("rs",data_plot_try$MarkerName))!=T)
# data_plot_try[jdx,1] <- paste0("chr",data_plot_try$CHR,":",data_plot_try$position)[jdx]
# kdx <- which(data_plot_try$CHR==discovery_snp$CHR.x[k]&
#                data_plot_try$position== discovery_snp$Pos[k])
# data_plot_try[kdx,1] <- paste0("chr",discovery_snp$CHR.x[k],":",discovery_snp$Pos[k])
# write.table(data_plot_try,file = paste0("./discovery_SNP/result/locus_plot_data/data_plot_",discovery_snp$SNP.ONCO[k],".txt"),sep="\t",col.names=T,row.names=F,quote=F)


#CIMBA analysis data
cimba_result_all <- as.data.frame(fread("./data/brca1_bcac_tn_meta.txt",header = T))
colnames(cimba_result_all)[10] <- "P"
cimba_result = cimba_result_all %>% 
  filter(Freq1>=0.008&
           Freq1<=0.992&
           CHR!=23) %>% 
  select(MarkerName,CHR,
         position,P)

colnames(cimba_result) <- c("MarkerName",
                            "CHR",
                            "position",
                            "P.value")
#rs17215231
k <- 29
  idx.temp <- which(cimba_result$CHR==discovery_snp$CHR.x[k]&
                      cimba_result$position>= discovery_snp$Pos[k]-10^6&
                      cimba_result$position<= discovery_snp$Pos[k]+10^6)
  cimba_result_try <- cimba_result[idx.temp,]
  #if the SNP don't have rs_id
  #use the chr:position format
  jdx <- which(c(1:nrow(cimba_result_try)%in%grep("rs",cimba_result_try$MarkerName))!=T)
  cimba_result_try[jdx,1] <- paste0("chr",cimba_result_try$CHR,":",cimba_result_try$position)[jdx]
  write.table(cimba_result_try,file = paste0("./discovery_SNP/result/locus_plot_data/cimba_result_plot",discovery_snp$SNP.ONCO[k],".txt"),sep="\t",col.names=T,row.names=F,quote=F)
 #rs2464195
 
  k <- 8
  idx.temp <- which(cimba_result$CHR==discovery_snp$CHR.x[k]&
                      cimba_result$position>= discovery_snp$Pos[k]-10^6&
                      cimba_result$position<= discovery_snp$Pos[k]+10^6)
  cimba_result_try <- cimba_result[idx.temp,]
  #if the SNP don't have rs_id
  #use the chr:position format
  jdx <- which(c(1:nrow(cimba_result_try)%in%grep("rs",cimba_result_try$MarkerName))!=T)
  cimba_result_try[jdx,1] <- paste0("chr",cimba_result_try$CHR,":",cimba_result_try$position)[jdx]
  write.table(cimba_result_try,file = paste0("./discovery_SNP/result/locus_plot_data/cimba_result_plot",discovery_snp$SNP.ONCO[k],".txt"),sep="\t",col.names=T,row.names=F,quote=F)
  
  
  #mixed-effect model
meta_result_shared_1p <- as.data.frame(fread("/Users/haoyuzhang/GoogleDrive/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/meta_result_shared_1p_mixed.txt",header=T))
  mtop <- meta_result_shared_1p
  temp.str <- strsplit(mtop$rs_id,":")
  n <- nrow(mtop)
  rs_id <- rep("c",n)
  for(i in 1:n){
    rs_id[i] <- temp.str[[i]][1]
  }
  mtop$rs_id <- rs_id
  mtop.p <- mtop$p.value
  #just for easy coding, put mtop as ftop
  ftop = mtop
  ftop$subtypes.p <- ftop$p.value
  ftop =  ftop %>% 
    mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0),1E-20,subtypes.p))
  subtypes_gwas_result <- ftop %>%
    select(rs_id,CHR,position,subtypes.p.new)
  colnames(subtypes_gwas_result) <- c("MarkerName",
                                      "CHR",
                                      "position",
                                      "P.value")
  idx <- NULL
  for(k in 1:nrow(discovery_snp)){
    idx.temp <- which(subtypes_gwas_result$CHR==discovery_snp$CHR.x[k]&
                        subtypes_gwas_result$position>= discovery_snp$Pos[k]-10^6&
                        subtypes_gwas_result$position<= discovery_snp$Pos[k]+10^6)
    data_plot_try <- subtypes_gwas_result[idx.temp,]
    #if the SNP don't have rs_id
    #use the chr:position format
    jdx <- which(c(1:nrow(data_plot_try)%in%grep("rs",data_plot_try$MarkerName))!=T)
    data_plot_try[jdx,1] <- paste0("chr",data_plot_try$CHR,":",data_plot_try$position)[jdx]
    write.table(data_plot_try,file = paste0("./discovery_SNP/result/locus_plot_data/mtop_",discovery_snp$SNP.ONCO[k],".txt"),sep="\t",col.names=T,row.names=F,quote=F)
    
  }

  #fixed-effect model
  meta_result_shared_1p <- as.data.frame(fread("/Users/haoyuzhang/GoogleDrive/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/meta_result_shared_1p_fixed.txt",header=T))
  ftop <- meta_result_shared_1p
  ftop.p <- ftop$p.value
  subtypes.p <- ftop.p
  ftop$subtypes.p <- ftop$p.value
  temp.str <- strsplit(ftop$rs_id,":")
  n <- nrow(ftop)
  rs_id <- rep("c",n)
  for(i in 1:n){
    rs_id[i] <- temp.str[[i]][1]
  }
  ftop$rs_id <- rs_id
  ftop =  ftop %>% 
    mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0),1E-20,subtypes.p))
  subtypes_gwas_result <- ftop %>%
    select(rs_id,CHR,position,subtypes.p.new)
  colnames(subtypes_gwas_result) <- c("MarkerName",
                                      "CHR",
                                      "position",
                                      "P.value")
  idx <- NULL
  for(k in 1:nrow(discovery_snp)){
    idx.temp <- which(subtypes_gwas_result$CHR==discovery_snp$CHR.x[k]&
                        subtypes_gwas_result$position>= discovery_snp$Pos[k]-10^6&
                        subtypes_gwas_result$position<= discovery_snp$Pos[k]+10^6)
    data_plot_try <- subtypes_gwas_result[idx.temp,]
    #if the SNP don't have rs_id
    #use the chr:position format
    jdx <- which(c(1:nrow(data_plot_try)%in%grep("rs",data_plot_try$MarkerName))!=T)
    data_plot_try[jdx,1] <- paste0("chr",data_plot_try$CHR,":",data_plot_try$position)[jdx]
    write.table(data_plot_try,file = paste0("./discovery_SNP/result/locus_plot_data/ftop_",discovery_snp$SNP.ONCO[k],".txt"),sep="\t",col.names=T,row.names=F,quote=F)
    
  }
  