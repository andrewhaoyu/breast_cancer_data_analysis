#purpose debug
#figure out where is wrong
#why the prs continues increase as more snps get in



setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction')

#load("./EB_whole_genome/result/whole_gonome.rdata")
load("./intrinsic_subtypes_whole_genome/ICOG/result/whole_gonome.rdata")

#######prs file need A1 to be effect_allele
#######we need to reverse all the log odds ratio since we put A2 as effect allele
#whole_genome = whole_genome %>% mutate(p.min = pmin(p.value,FTOP_result))
#head(whole_genome)
#save(whole_genome,file = "./EB_whole_genome/result/whole_gonome.rdata")
#LD clumping based on the min-p value of two-stage model and standard analysis
library(data.table)

library(dplyr)
#check duplicated
known_snp <- read.csv("/data/zhangh24/breast_cancer_data_analysis/data/212_known_discovery_snp_paper_order.csv")
known_snp = known_snp %>% 
  mutate(chr.pos = paste0(CHR,":",position))
whole_genome = whole_genome %>% 
  mutate(chr.pos = paste0(CHR,":",position))


idx <- which(whole_genome$chr.pos=="17:7571752")
whole_genome[idx,]
known_snp_all <- left_join(known_snp,whole_genome,
                                by="chr.pos")

known_snp_all[c(61,62,119,120),]
#take out duplicated snp
known_snp_all = known_snp_all[-c(61,120),]
idx <- which(duplicated(known_snp_all$Best.published.SNP))
length(idx)
#save(whole_genome_clump,file = "./EB_whole_genome/result/whole_genome_clump.rdata")
#No need to rerun the previous code again

#create prs files based on different p-threshold
setwd('/data/zhangh24/breast_cancer_data_analysis/risk_prediction')
library(dplyr)
library(data.table)
#load("./EB_whole_genome/result/whole_genome_clump.rdata")
dim(known_snp_all)
#method <- c("standard","two-stage","eb")
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,
            1E-04,5E-04,1E-03,5E-03,1E-02)

#create the prs file for two-stage and eb
subtypes <- c("Luminal_A",
              "Luminal_B",
              "Luminal_B_HER2Neg",
              "HER2_Enriched",
              "TN")
select.names <- subtypes
score <- known_snp_all %>%  select(select.names)

known_snp_new = known_snp_all %>% mutate(SNP=SNP.ONCO) %>% 
  select(Best.published.SNP,SNP,effect_allele,p.min) %>% 
  cbind(score)
known_snp_new = known_snp_new[complete.cases(known_snp_new),]

#generate extract snps list

extract_snp = known_snp_new %>% select(SNP)

write.table(extract_snp,file = paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extrac_snp_list.txt"),row.names=F,col.names=T,quote=F)



#generate map files 
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome/ONCO/ERPRHER2GRADE_fixed_baseline/result/onco_info.Rdata")
#map file have four columns, no header
#chromosome (1-22, X, Y or 0 if unplaced)
#rs# or snp identifier
#Genetic distance (morgans)
#Base-pair position (bp units)
onco_info$genetic_distance <- 0
map_file <- onco_info %>% 
  select(CHR,rs_id,genetic_distance,position)
write.table(map_file,
            file = "/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/map.txt",
            row.names = F,
            col.names=F,
            quote = F)

extract_rs_id <- rep("c",nrow(prs))
for(i in 1:nrow(prs)){
  extract_rs_id[i] <- as.character(prs[i,1])
}
temp <- which(extract_rs_id%in%
                as.character(map_file$rs_id))



prs.code <- rep("c",length(select.names))
prs.code <- data.frame(prs.code,stringsAsFactors=F)
temp <- 1




rm(list=ls())
#args=(commandArgs(TRUE))
#for(p in 1:length(args)){
#       eval(parse(text=args[[p]]))
#  }
#print(i)
#i1 <- i
#commandarg <- commandArgs(trailingOnly=F)
#myarg <- commandarg[length(commandarg)]
#myarg <- sub("-","",myarg)
#i <- as.numeric(myarg)
#print(i)
#pheno is ICOGS,data2 is onco_array
arg <- commandArgs(trailingOnly=T)
i <- as.numeric(arg[[1]])
i1 <- i
print(i)
library(R.utils)
setwd("/data/zhangh24/breast_cancer_data_analysis/")

subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))
pheno.file <- "./data/pheno.onco"
load(pheno.file)
n.sub = nrow(pheno)
y.pheno.mis2 <- cbind(pheno$Behaviour1,pheno$PR_status1,pheno$ER_status1,pheno$HER2_status1)
colnames(y.pheno.mis2) = c("Behaviour1","PR_status1",
                           "ER_status1","HER2_status1")


idx.fil <- onco.order[,1]%in%pheno$Onc_ID
n <- length(idx.fil)
snpvalue <- rep(0,n)
idx.match <- match(pheno$Onc_ID,onco.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
#load("./whole_genome/ONCO/ERPRHER2_fixed/result/score.test.support.onco.ERPRHER2.Rdata")




Filesdir <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_imputed"
Files <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.",full.names=T)
Filesex <- dir(Filesdir,pattern="OncoArray_european_merged_b1_15.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]


qctool.command <- rep("c",567)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)


for(i in 1:567){
  geno.file <- Files[i]
  
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g ",Files[i]," -incl-rsids /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extrac_snp_list.txt -og /data/zhangh24/BCAC/impute_onco/extracted_snp_test_",i,".txt")
  qctool.command[i,1] <- temp
  
}


write.table(qctool.command,file = paste0("/data/zhangh24/BCAC/impute_onco/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)





setwd("/data/zhangh24/breast_cancer_data_analysis/")


n.raw <- 142273
snpvalue <- rep(0,n.raw) 
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/onco_order.txt.gz"
onco.order <- read.table(gzfile(subject.file))


library(data.table)


load(paste0("./risk_prediction/result/split.id.rdata"))
#icog.test.id <- Generatetestid(subtypes.icog)
#load the sample data
icog.train.id <- split.id[[1]]
onco.train.id <- split.id[[2]]
onco.test.id <- split.id[[3]]
icog.vad.id <- split.id[[4]]
onco.vad.id <- split.id[[4]]
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")











#take out the first row of the sample data
#the first row of the sample data is the dataformat of each column
sample.data <- read.table("/data/zhangh24/BCAC/impute_onco/sample.txt",header=T,stringsAsFactors = F)
sample.data <- sample.data[-1,,drop=F]
n <- nrow(sample.data)  
extract.num <- nrow(known_snp_new)
snpid.result <- rep("c",extract.num)
n.sub <- nrow(sample.data)
snpvalue.result <- matrix(0,n.sub,extract.num)


total <- 0

for(i in 1:567){
  
  print(i)  
  geno.file <- paste0("/data/zhangh24/BCAC/impute_onco/extracted_snp_test_",i,".txt")
  num <- as.numeric(system(paste0('cat ',geno.file,' | wc -l '),intern=T))
  
  if(num!=0){
    
    
    con <- file(geno.file)
    temp <- 0
    open(con)
    for(i in 1:num){
      
      oneLine <- readLines(con,n=1)
      myVector <- strsplit(oneLine," ")
      snpid <- as.character(myVector[[1]][3])
      
      
      temp <- temp+1
      snpid.result[temp+total] <- snpid
      snpvalue <- rep(0,n)
      
      
      snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
      
      snpvalue <- 2-convert(snppro,n.raw)
      #snpvalue <- snpvalue[idx.fil][idx.match]
      snpvalue.result[,temp+total] <- snpvalue
      
      
    }
    close(con)
    
    total <- total+num
    
    
  }
  
  
}

snpid.result <- snpid.result[1:total]
snpvalue.result <- snpvalue.result[,1:total]
# snpvalue.result.order= snpvalue.result[idx.fil,]
#extract.list.shared <- prs[,1]
#temp2 = x.snp.j.order[,7]
#temp1 = snpvalue.result.order[,189]
#head(temp1)
#head(temp2)
#sum(abs(temp2-temp1))
#log.odds.intrinsic.all[189,]
#known_snp_new[202,]
#colnames(x.snp.j)[189]
extract.list.shared <- known_snp_new[complete.cases(known_snp_new[,2]),2]

idx.match <- match(extract.list.shared,snpid.result)
snpid.result <- snpid.result[idx.match]
all.equal(snpid.result,extract.list.shared)
snpvalue.result <- snpvalue.result[,idx.match]

colnames(snpvalue.result) <- snpid.result
n.snp <- 207
#prs.la <-   snpvalue.result[,1:n.snp]%*%prs[1:n.snp,3]/((n.snp+1)*2)
prs.la <-   (2-snpvalue.result[,1:n.snp])%*%known_snp_new[1:n.snp,5]/(n.snp*2)
#prs.la <-   snpvalue.result[,n.snp]*prs[n.snp,3]/(2)
prs.tn <- snpvalue.result%*%prs[1:n.snp,4]/(n.snp*2)
i <- 8
j <- 1
prs.plink <- read.table(paste0("/data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/test_out.profile"),header=T)
all.equal(as.numeric(prs.la),as.numeric(prs.plink[,4]))


load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/meta_result_shared_1p_082119.Rdata")

write.csv(prs,file= "/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/tempprs.csv")

colnames(meta_result_shared_1p)[c(16:20)] <- paste0("beta_",subtypes)
meta_result_shared_1p_sub = meta_result_shared_1p %>% 
  select(SNP.ONCO,beta_Luminal_A,beta_TN)
colnames(meta_result_shared_1p_sub)[1] <- "SNP"
prs_compare = left_join(prs,meta_result_shared_1p_sub,
                        by="SNP")
prs.la.temp <-   snpvalue.result[,1:n.snp]%*%prs_compare[1:n.snp,5]
prs.plink.temp <- prs.plink
prs.plink[,4] <- prs.la


# extract.result <- list(snpid.result,snpvalue.result)
# save(extract.result,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_result_shared.Rdata")
# 
# extract.list.shared <- extract.list[1:total,]
# idx.match <- match(extract.list.shared$SNP.ONCO,snpid.result)
# snpid.result <- snpid.result[idx.match]
# all.equal(snpid.result,extract.list.shared$SNP.ONCO)
# snpvalue.result <- snpvalue.result[,idx.match]
# extract.result <- list(snpid.result,snpvalue.result)
# colnames(snpvalue.result) <- snpid.result
# 
# 
# 
# extract.result <- list(snpid.result,snpvalue.result)
# save(extract.result,file="/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/ERPRHER2GRADE_fixed_baseline/result/extract_result_shared.Rdata")
# 
# 





prs.code <- paste0("/data/zhangh24/plink --dosage /data/zhangh24/BCAC/impute_onco_dosage/dosage_all noheader skip0=1 skip1=1 format=1 --map  /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/map.txt --fam /data/zhangh24/BCAC/impute_onco/onco_plink.fam --extract /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extrac_snp_list --out /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extract_snps_genotype")


"/data/zhangh24/qctool_v1.4-linux-x86_64/qctool -g chr22_test.gen -incl-rsids /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extrac_snp_list.txt -og /data/zhangh24/BCAC/impute_onco/extracted_snp_test"


"/data/zhangh24/software/qctool/qctool -g /data/zhangh24/BCAC/impute_onco/onco_all.gen.gz -incl-rsids /data/zhangh24/breast_cancer_data_analysis/risk_prediction/subtypes_prs/result/extrac_snp_list.txt -og /data/zhangh24/BCAC/impute_onco/extracted_snp_test"   
/spin1/users/zhangh24/software/qctool/qctool


