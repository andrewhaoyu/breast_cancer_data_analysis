#install_github("andrewhaoyu/bc2",ref='development', args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
#install_github("andrewhaoyu/bc2", args = c('--library="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.4"'))
###1 represent Icog
###2 represent Onco

rm(list=ls())
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i1 <- as.numeric(myarg)
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

library(readr)
library(devtools)
library(CompQuadForm)
library(bc2)
library(data.table)
# load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ICOG/ERPRHER2_fixed/result/extract_list.Rdata")
# load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2_fixed/result/extract_result.Rdata")

#idx <- which(extract.result[[1]]=="rs372562666:1:120561314:G:A")


dim(extract.result[[2]])


discovery.snp.onco <- as.data.frame(fread("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis_check_LD/discovery_onco_data.csv",header=T))
#x.test.all.mis1 <- discovery.snp.icog


discovery_snp <- read.csv("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_snp_summary_new.csv",header=T)

data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
data2 <- as.data.frame(data2)
LDandDp <- function(snp1,snp2){
  LD <- cor(snp1,snp2)^2
  p1 <- mean(snp1)/2
  p2 <- mean(snp2)/2
  #p12 <- cov(snp1,snp2)+p1*p2
  d <- mean(snp1*snp2)/4-p1*p2
  d <- 2*d
  if(d<0){
    dmax <- min(p1*p2,(1-p1)*(1-p2))
  }else{
    dmax <- min(p1*(1-p2),(1-p1)*p2)
  }
  Dp <- d/dmax
  return(list(LD,Dp))
  
}


#######snp rs150157076:120586681:A:C and snp rs11249433
i.dis <- which(colnames(discovery.snp.onco)=="rs150157076:120586681:A:C")
i.known <- which(colnames(data2)=="rs2532263")
idx.control <- which(data2$Behaviour1==0)
ld.vec <- rep(0,315)
snp2 <- as.vector(data2[idx.control,i.known])
for(i in 1:315){
  snp1 <- discovery.snp.onco[idx.control,i]
  ld.vec[i] <- LDandDp(snp1,snp2)[[1]]
}
ld.vec <- ld.vec[-143]
ld.vec[67]
which.min(ld.vec)
mean(snp1)/2

LDandDp(snp1,snp2)













#######snp rs150157076:120586681:A:C and snp rs11249433
i.dis <- which(colnames(discovery.snp.onco)=="rs150157076:120586681:A:C")
i.known <- which(colnames(data2)=="rs11249433")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)





#######snp rs11264454:156153043:A:Gand snp rs11249433
i.dis <- which(colnames(discovery.snp.onco)=="rs11264454:156153043:A:G")
i.known <- which(colnames(data2)=="rs4971059")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
mean(snp2)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs11749176:44145931:T:A and snp rs10941679
i.dis <- which(colnames(discovery.snp.onco)=="rs11749176:44145931:T:A")
i.known <- which(colnames(data2)=="rs10941679")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
mean(snp2)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)



#######snp c5_pos45369617 and snp rs10941679
i.dis <- which(colnames(discovery.snp.onco)=="c5_pos45369617")
i.known <- which(colnames(data2)=="rs10941679")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
mean(snp2)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs149663829:68152587:C:A and snp rs35951924
i.dis <- which(colnames(discovery.snp.onco)=="rs149663829:68152587:C:A")
i.known <- which(colnames(data2)=="rs35951924")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2

snp2 <- as.vector(data2[idx.control,i.known])
mean(snp2)/2
LDandDp(snp1,snp2)



#######snp rs6860806:131640536:A:G and snp rs6596100
i.dis <- which(colnames(discovery.snp.onco)=="rs6860806:131640536:A:G")
i.known <- which(colnames(data2)=="rs6596100")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2

snp2 <- as.vector(data2[idx.control,i.known])
mean(snp2)/2
LDandDp(snp1,snp2)


#######snp rs7760611 and snp rs2223621
i.dis <- which(colnames(discovery.snp.onco)=="rs7760611")
i.known <- which(colnames(data2)=="rs2223621")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2

snp2 <- as.vector(data2[idx.control,i.known])
mean(snp2)/2
LDandDp(snp1,snp2)

chr8_116679547_A_G
#######snp chr8_116679547_A_G and snp rs13267382
i.dis <- which(colnames(discovery.snp.onco)=="chr8_116679547_A_G")
i.known <- which(colnames(data2)=="rs13267382")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2

snp2 <- as.vector(data2[idx.control,i.known])
mean(snp2)/2
LDandDp(snp1,snp2)



#######snp rs12765365:64848937:T:C and snp rs10995201
i.dis <- which(colnames(discovery.snp.onco)=="rs12765365:64848937:T:C")
i.known <- which(colnames(data2)=="rs10995201")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2

snp2 <- as.vector(data2[idx.control,i.known])
mean(snp2)/2
LDandDp(snp1,snp2)



#######snp 12:29140260:G:A and snp rs7297051
i.dis <- which(colnames(discovery.snp.onco)=="12:29140260:G:A")
i.known <- which(colnames(data2)=="rs7297051")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2

snp2 <- as.vector(data2[idx.control,i.known])
mean(snp2)/2
LDandDp(snp1,snp2)


#######snp rs1061657 and snp rs1292011
i.dis <- which(colnames(discovery.snp.onco)=="rs1061657")
i.known <- which(colnames(data2)=="rs1292011")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2

snp2 <- as.vector(data2[idx.control,i.known])
mean(snp2)/2
LDandDp(snp1,snp2)

#######snp 17:43681771:C:T and snp rs2532263
i.dis <- which(colnames(discovery.snp.onco)=="17:43681771:C:T")
i.known <- which(colnames(data2)=="rs2532263")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2

snp2 <- as.vector(data2[idx.control,i.known])
mean(snp2)/2
LDandDp(snp1,snp2)
#17:43681771:C:T




#######snp rs6697258:120485335:C:A and snp rs11249433
i.dis <- which(colnames(discovery.snp.onco)=="rs6697258:120485335:C:A")
i.known <- which(colnames(data2)=="rs11249433")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)

#######snp 1:145126177:G:A and snp rs12405132
i.dis <- which(colnames(discovery.snp.onco)=="1:145126177:G:A")
i.known <- which(colnames(data2)=="rs12405132")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs348196:155666961:T:C and snp rs4971059
i.dis <- which(colnames(discovery.snp.onco)=="rs348196:155666961:T:C")
i.known <- which(colnames(data2)=="rs4971059")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)



#######snp rs11749176:44145931:T:A and snp rs10941679
i.dis <- which(colnames(discovery.snp.onco)=="rs11749176:44145931:T:A")
i.known <- which(colnames(data2)=="rs10941679")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)



#######snp 1:145126177 and snp rs12405132
i.dis <- which(colnames(discovery.snp.onco)=="rs56826596:45374890:G:A")
i.known <- which(colnames(data2)=="rs10941679")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)




#######snp 1:145126177 and snp rs12405132
i.dis <- which(colnames(discovery.snp.onco)=="1:145126177:G:A")
i.known <- which(colnames(data2)=="rs12405132")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
mean(snp1)/2
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs6677545 and snp rs35383942
i.dis <- which(colnames(discovery.snp.onco)=="rs6677545:200342046:A:C")
i.known <- which(colnames(data2)=="rs35383942")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)

#######snp rs6677545 and snp rs6678914
i.dis <- which(colnames(discovery.snp.onco)=="rs6677545:200342046:A:C")
i.known <- which(colnames(data2)=="rs6678914")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs56826596 and snp rs10941679
i.dis <- which(colnames(discovery.snp.onco)=="rs56826596:45374890:G:A")
i.known <- which(colnames(data2)=="rs10941679")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)

#######snp rs139331653:45939294:G:A and snp rs10941679
i.dis <- which(colnames(discovery.snp.onco)=="rs139331653:45939294:G:A")
i.known <- which(colnames(data2)=="rs10941679")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)

#######snp rs34044188:65257363:C:T and snp rs10995201
i.dis <- which(colnames(discovery.snp.onco)=="rs34044188:65257363:C:T")
i.known <- which(colnames(data2)=="rs10995201")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs17743054:42900892:T:C and snp rs6507583
i.dis <- which(colnames(discovery.snp.onco)=="rs17743054:42900892:T:C")
i.known <- which(colnames(data2)=="rs6507583")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)


#######snp rs16988381 and snp rs17879961
i.dis <- which(colnames(discovery.snp.onco)=="rs16988381")
i.known <- which(colnames(data2)=="rs17879961")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)

#######snp rs16988381 and snp rs132390
i.dis <- which(colnames(discovery.snp.onco)=="rs16988381")
i.known <- which(colnames(data2)=="rs132390")
idx.control <- which(data2$Behaviour1==0)
snp1 <- discovery.snp.onco[idx.control,i.dis]
snp2 <- as.vector(data2[idx.control,i.known])
LDandDp(snp1,snp2)

