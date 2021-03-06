#test mahattan plot package to see which one looks best
library(qqman)
manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P",ylim =c(0,12))
qq(gwasResults$P,ylim  = c(0,12))
library(data.table)
library(dplyr)
data <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))

gwas_result <- data %>%
  #filter(p.meta<=10^-4) %>% 
  select(c(SNP.Onco,chr.Onco,Position.Onco,p.meta))
colnames(gwas_result) <- c("SNP",
                           "CHR",
                           "BP",
                           "P")
Calculatelambda <- function(x,stat_type){
  if(stat_type!="CHISQ4"&
     stat_type!="CHISQ5"){
    if (stat_type == "Z")
      z = x
    
    if (stat_type == "CHISQ")
      z = sqrt(x)
    
    if (stat_type == "PVAL")
      z = qnorm(x / 2)
    lambda = round(median(z^2) / qchisq(0.5,1), 3)
    return(lambda)
    
  }else{
    if(stat_type =="CHISQ5"){ z = qchisq(x,5,lower.tail = F)
    lambda = round(median(z) / qchisq(0.5,5), 3)}
    if(stat_type =="CHISQ4"){  z = qchisq(x,4,lower.tail = F)
    lambda = round(median(z) / qchisq(0.5,4), 3)}
    return(lambda)
  }
  
  }

Calculatelambda(gwas_result$P,"PVAL")
Calculatelambda(mtop.p,"CHISQ4")
Calculatelambda(ftop.p,"CHISQ5")
gwas_result = gwas_result %>%
 mutate(P= ifelse(P==0,1E-20,P))


fine_mapping <- read.csv("./data/filter_regions_standard.csv",header= T)
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
dim(gwas_result_filter)
min(gwas_result_filter$P)

gwas_result_filter[idx,]
idx.sig <- which(gwas_result_filter$P<=5*10^-8)
length(idx.sig)
head(gwas_result_filter[idx.sig,])
gwas_sig <- gwas_result_filter[idx.sig,]
head(gwas_result_filter[idx.sig,])
gwas_sig[150:550,]
manhattan(gwas_result, col= c("blue4", "orange3"),
          cex = 0.6,suggestiveline = F)
qq(gwas_result$P)
library(GenABEL)

estlambda(Pvalues, method="median")
library("CMplot")
CMplot(gwas_result_filter,plot.type="c")

#######mahattan plot for subtypes analysis
load("./data/whole_genome_ftop_1p.Rdata")
ftop <- meta_result_shared_1p
load("./data/whole_genome_mtop_1p.Rdata")
mtop <- meta_result_shared_1p
ftop.p <- ftop$p.value
ftop$subtypes.p <- ftop$p.value
ftop =  ftop %>% 
  mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0),1E-20,subtypes.p))
subtypes_gwas_result <- ftop %>%
  filter(subtypes.p<=1E-04) %>% 
  select(SNP.ONCO,CHR,position,subtypes.p.new)

colnames(subtypes_gwas_result) <- c("SNP",
                                    "CHR",
                                    "BP",
                                    "P")
manhattan(subtypes_gwas_result, col= c("blue4", "orange3"),
          cex = 0.6,suggestiveline = F)  


mtop$subtypes.p <- mtop$p.value
mtop =  mtop %>% 
  mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0),1E-20,subtypes.p))
subtypes_gwas_result <- ftop %>%
  filter(subtypes.p<=1E-04) %>% 
  select(SNP.ONCO,CHR,position,subtypes.p.new)

colnames(subtypes_gwas_result) <- c("SNP",
                                    "CHR",
                                    "BP",
                                    "P")
manhattan(subtypes_gwas_result, col= c("blue4", "orange3"),
          cex = 0.6,suggestiveline = F) 








subtypes.p <- apply(cbind(ftop.p,mtop.p),1,min)
subtypes.p <- ifelse(is.na(subtypes.p),1E-20,subtypes.p)
ftop.p <- ifelse(is.nan(ftop.p),1E-20,ftop.p)
Calculatelambda(ftop.p,"PVAL")
Calculatelambda(mtop.p,"PVAL")


idx <- which(is.nan(ftop$subtypes.p))
ftop[idx,]


ftop$subtypes.p <- subtypes.p
ftop =  ftop %>% 
  mutate(subtypes.p.new = ifelse((is.nan(subtypes.p)==T)|(subtypes.p==0),1E-20,subtypes.p))
subtypes_gwas_result <- ftop %>%
  filter(subtypes.p<=1E-04) %>% 
  select(SNP.ONCO,CHR,position,subtypes.p.new)

colnames(subtypes_gwas_result) <- c("SNP",
                                           "CHR",
                                           "BP",
                                           "P")

fine_mapping <- read.csv("./data/filter_regions_subtypes.csv",header= T)
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
  idx <- which(subtypes_gwas_result$CHR==chr_temp&subtypes_gwas_result$BP>=start_temp&subtypes_gwas_result$BP<=end_temp)
  idx_cut <- c(idx_cut,idx)
}
idx_cut <- unique(idx_cut)
subtypes_gwas_result_filter <- subtypes_gwas_result[-idx_cut,]
manhattan(subtypes_gwas_result_filter, col= c("blue4", "orange3"),
          cex = 0.6,suggestiveline = F)     
manhattan(subtypes_gwas_result, col= c("blue4", "orange3"),
          cex = 0.6,suggestiveline = F)   
gwas_sig <- subtypes_gwas_result_filter %>% filter(P<=5E-08)
gwas_sig[order(gwas_sig$CHR),]


#####mahattan plot for CIMBA BRCA1 META
library(data.table)
# cimba_result <- fread("./data/brca1_bcac_tn_ma1.txt",header = T) 

# write.table(cimba_result,file = "./data/brca1_bcac_tn_meta.txt",row.names = F,
#             col.names=T,quote=F)
cimba_result_all <- as.data.frame(fread("./data/brca1_bcac_tn_meta.txt",header = T))
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
  filter(P<=1E-04&
           Freq1>=0.008&
           Freq1<=0.992) %>% 
  select(MarkerName,CHR,
         position,P)

colnames(cimba_result) <- c("SNP",
                            "CHR",
                            "BP",
                            "P")
fine_mapping <- read.csv("./data/filter_regions_cimba.csv",header= T)
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
  idx <- which(cimba_result$CHR==chr_temp&cimba_result$BP>=start_temp&cimba_result$BP<=end_temp)
  idx_cut <- c(idx_cut,idx)
}
idx_cut <- unique(idx_cut)
cimba_result_filter <- cimba_result[-idx_cut,]
manhattan(cimba_result_filter, col= c("blue4", "orange3"),
          cex = 0.6,suggestiveline = F)     

gwas_sig <- cimba_result_filter %>% filter(P<=5E-08)
gwas_sig[order(gwas_sig$CHR),]





gwas_result_filter$chr.pos <- paste0(gwas_result_filter$CHR,":",gwas_result_filter$BP)
subtypes_gwas_result_filter$chr.pos <- paste0(subtypes_gwas_result_filter$CHR,":",subtypes_gwas_result_filter$BP)
cimba_result_filter$chr.pos <- paste0(cimba_result_filter$CHR,":",cimba_result_filter$BP)

new.data <- merge(gwas_result_filter,
                  subtypes_gwas_result_filter,
                  by.x = "chr.pos",
                  by.y = "chr.pos",
                  all.x = T,
                  all.y = T)
new.data2 <- merge(new.data,
                   cimba_result_filter,
                   by.x="chr.pos",
                   by.y = "chr.pos",
                   all.x = T,
                   all.y = T)
n.snp <- nrow(new.data2)
temp <- strsplit(new.data2$chr.pos,":")
CHR <- rep(0,n.snp)
BP <- rep(0,n.snp)
for(i in 1:n.snp){
  CHR[i] <- as.numeric(temp[[i]][1])
  BP[i] <- as.numeric(temp[[i]][2])
}
new.data2$CHR <- CHR
new.data2$BP <- BP

plot.data = new.data2 %>% 
  select(chr.pos,CHR,BP,P.x,P.y,P)
colnames(plot.data) <- c("SNP",
                         "CHR","BP","Trait1",
                         "Trait2","Trait3")


library(CMplot)
data(pig60K)
setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/discovery_SNP/result")
CMplot(plot.data, plot.type="c", chr.labels=paste("Chr",c(1:22,"X"),sep=""), r=0.4, cir.legend=TRUE,
       outward=FALSE, cir.legend.col="black", cir.chr.h=1.3 ,chr.den.col="black", file="jpg",
       memo="", dpi=300)
