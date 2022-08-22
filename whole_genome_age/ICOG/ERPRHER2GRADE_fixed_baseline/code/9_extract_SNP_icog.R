
library(R.utils)
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")

n <- 109713
snpvalue <- rep(0,n)
subject.file <- "/data/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))

setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
x.covar.mis1 <- data1[,c(5:14,204)]
age <- data1[,204]
idx.complete <- which(age!=888)
y.pheno.mis1 <- y.pheno.mis1[idx.complete,]
x.covar.mis1 <- x.covar.mis1[idx.complete,]
SG_ID <- SG_ID[idx.complete]

idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
#library(bc2)
#load("./whole_genome_age/ICOG/ERPRHER2_fixed/result/score.test.support.icog.ERPRHER2.Rdata")

######### you need to modify this section#############################
#########load data_clean_select_sig in the investigate results code
load("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_list.Rdata")
#########snp.icogs.extract.id is the column SNP.ICOGS
#find the corresponding column for SNP.ICOG in data_clean_select_sig
#your outcome format is different with the old version, you can just do snp.icogs.extract.id = data_clean_select_sig$SNP.ICOGS
snp.icogs.extract.id <- extract.list[,(ncol(extract.list)-2),drop=F] 
########change the directory to your local directory
write.table(snp.icogs.extract.id,file = paste0("./whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_id_icog.txt"),quote = F,row.names=F)
######################################################################

#we don't need CHR 23 during extrac steps since there is no significant SNPs
Filesdir <- "/data/NC_BW/icogs_onco/genotype/imputed2/icogs_imputed"
Files <- dir(Filesdir,pattern="icogs_merged_b1_12.",full.names=T)
Filesex <- dir(Filesdir,pattern="icogs_merged_b1_12.chr23",full.names=T)
idx.sex <- Files%in%Filesex
Files <- Files[!idx.sex]
library(gtools)
Files <- mixedsort(Files)



qctool.command <- rep("c",564)
qctool.command <- data.frame(qctool.command,stringsAsFactors=F)

#we use the software qctool to extract the targeted SNPs from all the genotype files
#qctool is a command line software
#more details of the software can be found here: https://www.well.ox.ac.uk/~gav/qctool_v1/#overview
#to download the software: https://www.well.ox.ac.uk/~gav/qctool_v1/#download
#use the linux-x86_64 version
#put the software to your local folder to use
#qc tool is a command line software
#-g put into the genotype data
#-incl-rsids put the target gentoype id
#-og put the outcome file
#change the direcotry to your local directory
for(i in 1:564){
  geno.file <- Files[i]
  temp <- paste0("/data/zhangh24/qctool_v1.4-linux-x86_64/qctool ",
                 "-g ",Files[i]," ",
                 "-incl-rsids /data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/extract_id_icog.txt ",
                 "-og /data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/result/ERPRHER2_extract",i,".txt")
  qctool.command[i,1] <- temp
  
}
#change the directory to your local directory
write.table(qctool.command,file = paste0("/data/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_fixed_baseline/code/qc_extract_snp.sh"),col.names = F,row.names = F,quote=F)
#write qctool.command to a .sh file to submit to biowulf for running
#you can do swarm -f qc_extract_snp.sh -g 20 --time 4:00:00

