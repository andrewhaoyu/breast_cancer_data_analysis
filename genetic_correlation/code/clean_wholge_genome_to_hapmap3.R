setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/")
library(data.table)
hap3 <- fread("./data/HM3snplist_withCHRBP.txt",header=T)
dim(hap3)

CHR_generate <- function(position){
  n <- length(position)
  CHR <- rep(22,n)
  position.diff <- diff(position)
  idx <- which(position.diff<0)
  start <- 1
  temp <- 1
  for(i in 1:length(idx)){
    
    end <- idx[i]
    CHR[start:end] <- temp
    start <- end+1
    temp = temp+1
  }
  return(CHR)
}

hap3.chr.bp <- cbind(hap3[,4],hap3[,5])

hap3.chr.bp.uni <- apply(hap3.chr.bp,1,paste0,collapse = ":")
dim(hap3)
length(hap3.chr.bp.uni)
hap3$chr.pos <- hap3.chr.bp.uni
##############case control
##############ICOG

load("./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/ICOG/icog_result_odds_sd.Rdata")

CHR <- CHR_generate(icog_result$position)
icog_result.chr <- CHR
icog.chr.bp <- cbind(icog_result.chr,icog_result$position)
icog.chr.pos <- apply(icog.chr.bp,1,paste0,collapse=":")
icog_result$chr.pos <- icog.chr.pos
save(icog.chr.pos,file="./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/ICOG/icog.chr.pos.Rdata")

icog_new <- merge(hap3,icog_result,by.x="chr.pos",by.y = "chr.pos")

shared.data.icog <- icog_new[,c(1,2,5,6,8,10,18,19)]
head(shared.data.icog)
n <- nrow(shared.data.icog)
A1 <- rep(0,nrow(shared.data.icog))
A2 <- rep(0,nrow(shared.data.icog))
temp <- strsplit(shared.data.icog$rs_id,":")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}

shared.data.icog$A1 <- A1
shared.data.icog$A2 <- A2

shared.icog.data.complete <- na.omit(shared.data.icog)
save(shared.icog.data.complete,file="./genetic_correlation/standard_analysis/result/standard_gwas_result_hapmap3/case_control/icog/shared.icog.data.complete.Rdata")
###############ONCO

load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/ONCO/onco_result_odds_sd.Rdata")




CHR <- CHR_generate(onco_result$position)
onco_result.chr <- CHR
onco.chr.bp <- cbind(onco_result.chr,onco_result$position)
onco.chr.pos <- apply(onco.chr.bp,1,paste0,collapse=":")
onco_result$chr.pos <- onco.chr.pos
save(onco.chr.pos,file="./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/onco/onco.chr.pos.Rdata")

onco_new <- merge(hap3,onco_result,by.x="chr.pos",by.y = "chr.pos")

shared.data.onco <- onco_new[,c(1,2,5,6,8,10,18,19)]
head(shared.data.onco)
n <- nrow(shared.data.onco)
A1 <- rep(0,nrow(shared.data.onco))
A2 <- rep(0,nrow(shared.data.onco))
temp <- strsplit(shared.data.onco$rs_id,":")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}

shared.data.onco$A1 <- A1
shared.data.onco$A2 <- A2




shared.onco.data.complete <- na.omit(shared.data.onco)
save(shared.onco.data.complete,file="./genetic_correlation/standard_analysis/result/standard_gwas_result_hapmap3/case_control/ONCO/shared.onco.data.complete.Rdata")

##############meta 


load("./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/meta/case_control_meta_result_final.Rdata")
meta.chr.bp <- cbind(case_control_meta_result_final$CHR,
                     case_control_meta_result_final$position)
meta.chr.pos <- apply(bcac.chr.bp,1,paste0,collapse=":")

case_control_meta_result_final$chr.pos <- meta.chr.pos
save(meta.chr.pos,file="./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/meta/meta.chr.pos.Rdata")
shared.data <- merge(hap3,case_control_meta_result_final,by.x="chr.pos",by.y="chr.pos")
shared.data.meta <- shared.data[,c(1,2,5,6,7,10,11)]

head(shared.data.meta)
n <- nrow(shared.data.meta)
A1 <- rep(0,nrow(shared.data.meta))
A2 <- rep(0,nrow(shared.data.meta))
temp <- strsplit(shared.data.meta$rs_id,":")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}

shared.data.meta$A1 <- A1
shared.data.meta$A2 <- A2

shared.data.meta.complete <- na.omit(shared.data.meta)

shared.data.meta.complete.new <- merge(shared.data.meta.complete,shared.icog.data.complete,by.x="chr.pos",by.y="chr.pos")
shared.data.meta.complete <- shared.data.meta.complete.new[,c(1,2,3,4,5,14,6,7,8,9)]
colnames(shared.data.meta.complete) <- colnames(shared.icog.data.complete)
save(shared.data.meta.complete,file="./genetic_correlation/standard_analysis/result/standard_gwas_result_hapmap3/case_control/meta/shared.data.meta.complete.Rdata")











##############ER+_control
##############ICOG

load("./genetic_correlation/standard_analysis/result/standard_gwas_result/ER+_control/ICOG/icog_result_odds_sd.Rdata")

CHR <- CHR_generate(icog_result$position)
icog_result.chr <- CHR
icog.chr.bp <- cbind(icog_result.chr,icog_result$position)
icog.chr.pos <- apply(icog.chr.bp,1,paste0,collapse=":")
icog_result$chr.pos <- icog.chr.pos
save(icog.chr.pos,file="./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/ICOG/icog.chr.pos.Rdata")

icog_new <- merge(hap3,icog_result,by.x="chr.pos",by.y = "chr.pos")

shared.data.icog <- icog_new[,c(1,2,5,6,8,10,18,19)]
head(shared.data.icog)
n <- nrow(shared.data.icog)
A1 <- rep(0,nrow(shared.data.icog))
A2 <- rep(0,nrow(shared.data.icog))
temp <- strsplit(shared.data.icog$rs_id,":")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}

shared.data.icog$A1 <- A1
shared.data.icog$A2 <- A2

shared.icog.data.complete <- na.omit(shared.data.icog)
save(shared.icog.data.complete,file="./genetic_correlation/standard_analysis/result/standard_gwas_result_hapmap3/ER+_control/icog/shared.icog.data.complete.Rdata")
###############ONCO

load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/standard_analysis/result/standard_gwas_result/ER+_control/ONCO/onco_result_odds_sd.Rdata")




CHR <- CHR_generate(onco_result$position)
onco_result.chr <- CHR
onco.chr.bp <- cbind(onco_result.chr,onco_result$position)
onco.chr.pos <- apply(onco.chr.bp,1,paste0,collapse=":")
onco_result$chr.pos <- onco.chr.pos
save(onco.chr.pos,file="./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/onco/onco.chr.pos.Rdata")

onco_new <- merge(hap3,onco_result,by.x="chr.pos",by.y = "chr.pos")

shared.data.onco <- onco_new[,c(1,2,5,6,8,10,18,19)]
head(shared.data.onco)
n <- nrow(shared.data.onco)
A1 <- rep(0,nrow(shared.data.onco))
A2 <- rep(0,nrow(shared.data.onco))
temp <- strsplit(shared.data.onco$rs_id,":")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}

shared.data.onco$A1 <- A1
shared.data.onco$A2 <- A2




shared.onco.data.complete <- na.omit(shared.data.onco)
save(shared.onco.data.complete,file="./genetic_correlation/standard_analysis/result/standard_gwas_result_hapmap3/ER+_control/ONCO/shared.onco.data.complete.Rdata")

##############meta 


load("./genetic_correlation/standard_analysis/result/standard_gwas_result/ER+_control/meta/ER+_control_meta_result_final.Rdata")
# meta.chr.bp <- cbind(case_control_meta_result_final$CHR,
#                      case_control_meta_result_final$position)
# meta.chr.pos <- apply(bcac.chr.bp,1,paste0,collapse=":")

ERPos_control_meta_result_final$chr.pos <- meta.chr.pos
#save(meta.chr.pos,file="./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/meta/meta.chr.pos.Rdata")
shared.data <- merge(hap3,ERPos_control_meta_result_final,by.x="chr.pos",by.y="chr.pos")
shared.data.meta <- shared.data[,c(1,2,5,6,7,10,11)]

head(shared.data.meta)
n <- nrow(shared.data.meta)
A1 <- rep(0,nrow(shared.data.meta))
A2 <- rep(0,nrow(shared.data.meta))
temp <- strsplit(shared.data.meta$rs_id,":")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}

shared.data.meta$A1 <- A1
shared.data.meta$A2 <- A2

shared.data.meta.complete <- na.omit(shared.data.meta)

shared.data.meta.complete.new <- merge(shared.data.meta.complete,shared.icog.data.complete,by.x="chr.pos",by.y="chr.pos")
shared.data.meta.complete <- shared.data.meta.complete.new[,c(1,2,3,4,5,14,6,7,8,9)]
colnames(shared.data.meta.complete) <- colnames(shared.icog.data.complete)
save(shared.data.meta.complete,file="./genetic_correlation/standard_analysis/result/standard_gwas_result_hapmap3/ER+_control/meta/shared.data.meta.complete.Rdata")











##############ER-_control
##############ICOG

load("./genetic_correlation/standard_analysis/result/standard_gwas_result/ER-_control/ICOG/icog_result_odds_sd.Rdata")

CHR <- CHR_generate(icog_result$position)
icog_result.chr <- CHR
icog.chr.bp <- cbind(icog_result.chr,icog_result$position)
icog.chr.pos <- apply(icog.chr.bp,1,paste0,collapse=":")
icog_result$chr.pos <- icog.chr.pos
save(icog.chr.pos,file="./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/ICOG/icog.chr.pos.Rdata")

icog_new <- merge(hap3,icog_result,by.x="chr.pos",by.y = "chr.pos")

shared.data.icog <- icog_new[,c(1,2,5,6,8,10,18,19)]
head(shared.data.icog)
n <- nrow(shared.data.icog)
A1 <- rep(0,nrow(shared.data.icog))
A2 <- rep(0,nrow(shared.data.icog))
temp <- strsplit(shared.data.icog$rs_id,":")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}

shared.data.icog$A1 <- A1
shared.data.icog$A2 <- A2

shared.icog.data.complete <- na.omit(shared.data.icog)
save(shared.icog.data.complete,file="./genetic_correlation/standard_analysis/result/standard_gwas_result_hapmap3/ER-_control/icog/shared.icog.data.complete.Rdata")
###############ONCO

load("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/genetic_correlation/standard_analysis/result/standard_gwas_result/ER-_control/ONCO/onco_result_odds_sd.Rdata")




CHR <- CHR_generate(onco_result$position)
onco_result.chr <- CHR
onco.chr.bp <- cbind(onco_result.chr,onco_result$position)
onco.chr.pos <- apply(onco.chr.bp,1,paste0,collapse=":")
onco_result$chr.pos <- onco.chr.pos
save(onco.chr.pos,file="./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/onco/onco.chr.pos.Rdata")

onco_new <- merge(hap3,onco_result,by.x="chr.pos",by.y = "chr.pos")

shared.data.onco <- onco_new[,c(1,2,5,6,8,10,18,19)]
head(shared.data.onco)
n <- nrow(shared.data.onco)
A1 <- rep(0,nrow(shared.data.onco))
A2 <- rep(0,nrow(shared.data.onco))
temp <- strsplit(shared.data.onco$rs_id,":")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}

shared.data.onco$A1 <- A1
shared.data.onco$A2 <- A2




shared.onco.data.complete <- na.omit(shared.data.onco)
save(shared.onco.data.complete,file="./genetic_correlation/standard_analysis/result/standard_gwas_result_hapmap3/ER-_control/ONCO/shared.onco.data.complete.Rdata")

##############meta 


load("./genetic_correlation/standard_analysis/result/standard_gwas_result/ER-_control/meta/ER-_control_meta_result_final.Rdata")
# meta.chr.bp <- cbind(case_control_meta_result_final$CHR,
#                      case_control_meta_result_final$position)
# meta.chr.pos <- apply(bcac.chr.bp,1,paste0,collapse=":")

ERNeg_control_meta_result_final$chr.pos <- meta.chr.pos
#save(meta.chr.pos,file="./genetic_correlation/standard_analysis/result/standard_gwas_result/case_control/meta/meta.chr.pos.Rdata")
shared.data <- merge(hap3,ERNeg_control_meta_result_final,by.x="chr.pos",by.y="chr.pos")
shared.data.meta <- shared.data[,c(1,2,5,6,7,10,11)]

head(shared.data.meta)
n <- nrow(shared.data.meta)
A1 <- rep(0,nrow(shared.data.meta))
A2 <- rep(0,nrow(shared.data.meta))
temp <- strsplit(shared.data.meta$rs_id,":")
for(i in 1:n){
  print(i)
  
  A1[i] <- temp[[i]][3]
  A2[i] <- temp[[i]][4]
}

shared.data.meta$A1 <- A1
shared.data.meta$A2 <- A2

shared.data.meta.complete <- na.omit(shared.data.meta)

shared.data.meta.complete.new <- merge(shared.data.meta.complete,shared.icog.data.complete,by.x="chr.pos",by.y="chr.pos")
shared.data.meta.complete <- shared.data.meta.complete.new[,c(1,2,3,4,5,14,6,7,8,9)]
colnames(shared.data.meta.complete) <- colnames(shared.icog.data.complete)
save(shared.data.meta.complete,file="./genetic_correlation/standard_analysis/result/standard_gwas_result_hapmap3/ER-_control/meta/shared.data.meta.complete.Rdata")
