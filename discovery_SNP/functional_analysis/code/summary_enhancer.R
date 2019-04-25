#Goal: summary enhancer results from Jona
setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis")
data <- read.csv("./discovery_SNP/functional_analysis/result/enhancer_count.csv",header=T)
data <- data[,1:16]

discover <- read.csv("./discovery_SNP/additive_model/result/additive_model_result.csv")
colnames(discover)[1] <- "SNP"


try <- merge(data,discover,by.x= "Rs.id",by.y = "SNP")

new.sig <- ifelse(try$Mixed.Model.global.heterogeneity.test.p.value..baseline.ER.fixed..1<=0.05/32,1,0)
table(new.sig,data$Opposite.direction)

try$CCV_in_anyswitch_enh <- anyswitch

write.csv(try,file = "./discovery_SNP/functional_analysis/result/enhancer_count_update.csv")

totalswitch <- data$OFF.PRIMED+data$ACTIVE.OFF+data$ACTIVE.OFF.PRIMED+data$ACTIVE.PRIMED

anyswitch <- ifelse(totalswitch>=1,1,0)

table(data$Opposite.direction,data$CCV_in_switched_enh)

table(data$heterogeneous.test.results,anyswitch,data$Opposite.direction)

table(data$heterogeneous.test.results,data$CCV_in_switched_enh,data$Opposite.direction)

table(data$CCV_in_switched_enh,anyswitch)

table(data$heterogeneous.test.results,anyswitch)
as.matrix(table(anyswitch,data$Opposite.direction))
chisq.test(as.matrix(table(anyswitch,data$heterogeneous.test.results)))




data <- read.csv("./discovery_SNP/functional_analysis/result/enhancer_update_new.csv",header=T)
table(data$heterogeneous.test.results,
      data$CCV_in_ACTIVEswitched_enh)

result <- matrix("c",22,1)



for(i in 1:22){
  result[i,1] <- paste0("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
}

write.table(result,file = "./multi_ethnic/code/download_kg_vcf.sh",row.names=F,quote=F,col.names = F)

