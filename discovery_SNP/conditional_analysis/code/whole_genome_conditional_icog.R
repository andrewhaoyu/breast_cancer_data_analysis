args = commandArgs(trailingOnly = T)
i1 <- as.numeric(args[[1]])


load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")
library(data.table)
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")

n.raw <- 109713
snpvalue <- rep(0,n.raw)
subject.file <- "/gpfs/gsfs4/users/NC_BW/icogs_onco/genotype/imputed2/icogs_order.txt.gz"
Icog.order <- read.table(gzfile(subject.file))
setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data1 <- as.data.frame(data1)
y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
#x.test.all.mis1 <- data1[,c(27:206)]
SG_ID <- data1$SG_ID
x.covar.mis1 <- data1[,c(5:14,204)]
idx.fil <- Icog.order[,1]%in%SG_ID
idx.match <- match(SG_ID,Icog.order[idx.fil,1])
#Icog.order.match <- Icog.order[idx.fil,1][idx.match]
library(bc2)









text <- system(paste0("cat | wc -l /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_icog",i1,".txt"),intern=T)
extract.num <- as.integer(gsub(paste0(" /spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_icog",i1,".txt"),"",text))

if(extract.num==0){
  conditional.analysis.icog <- NULL
  save(conditional.analysis.icog,paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.analysis.icog",i1,".Rdata"))
  
}else{
  snpid.result <- rep("c",extract.num)
  score.result.icog <- 
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/all.conditional.snps.Rdata")
  setwd("/spin1/users/zhangh24/breast_cancer_data_analysis/")
  library(readr)
  library(devtools)
  library(CompQuadForm)
  library(bc2)
  library(data.table)
  library(bigmemory)
  load("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/support.matrix.Rdata")
  z.standard <- support.matrix[[1]]
  z.additive.design <- support.matrix[[2]]
  M <- support.matrix[[3]]
  number.of.tumor <- support.matrix[[4]]
  z.design.score.baseline <- support.matrix[[5]]
  z.design.score.casecase <- support.matrix[[6]]
  z.design.score.baseline.ER <- support.matrix[[7]]
  z.design.score.casecase.ER <- support.matrix[[8]]
  n.condition <- nrow(all.conditional.snps)
  
  data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
  data1 <- as.data.frame(data1)
  age1 <- data1[,204]
  idx.complete1 <- which(age1!=888)
  age1 <- age1[idx.complete1]
  
  y.pheno.mis1 <- cbind(data1$Behaviour1,data1$ER_status1,data1$PR_status1,data1$HER2_status1,data1$Grade1)
  colnames(y.pheno.mis1) = c("Behavior","ER","PR","HER2","Grade")
  
  known.all.mis1 <- data1[,c(27:203)]
  ###fake the onco array only snp
  
  
  
  x.covar.mis1 <- data1[,c(5:14)]
  
  icog.julie <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_icog.csv")
  icog.julie <- icog.julie[,-1]
  discovery.snp.icog <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_icog_data.csv",header=T)
  onco.julie <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/Julie_snp_onco.csv")
  onco.julie <- onco.julie[,-1]
  discovery.snp.onco <- fread("/spin1/users/zhangh24/breast_cancer_data_analysis/data/discovery_onco_data.csv")
  x.discovery.mis1 <- as.data.frame(cbind(icog.julie,discovery.snp.icog))
  x.discovery.mis2 <- as.data.frame(cbind(onco.julie,discovery.snp.onco))
  ###two snps onco array only
  sudo.icog.na <- rep(NA,nrow(data1))
  
  known.all.mis1 <- cbind(known.all.mis1,sudo.icog.na,x.discovery.mis1)
  known.all.mis2 <- cbind(known.all.mis2,x.discovery.mis2)
  
  y.pheno.mis1 <- y.pheno.mis1[idx.complete1,]
  x.covar.mis1 <- x.covar.mis1[idx.complete1,]
  x.covar.mis1 <- cbind(x.covar.mis1,age1)
  
  known.all.mis1 <- known.all.mis1[idx.complete1,]
  y.pheno.mis2 <- y.pheno.mis2[idx.complete2,]
  x.covar.mis2 <- x.covar.mis2[idx.complete2,]
  x.covar.mis2 <- cbind(x.covar.mis2,age2)
  
  known.all.mis2 <- known.all.mis2[idx.complete2,]
  
  known.flag.all <- all.conditional.snps$known.flag
  
  
  total <- 0
  
  for(i in i1:i1){
    
    
    print(i)  
    geno.file <- paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional_extract_icog",i,".txt"
    )
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
        #snpvalue <- rep(0,n)
        
        
        snppro <- as.numeric(unlist(myVector)[7:length(myVector[[1]])])
        
        snpvalue <- convert(snppro,n.raw)
        snpvalue <- snpvalue[idx.fil][idx.match]
        snpvalue.result[,temp+total] <- snpvalue
        
        
      }
      close(con)
      
      total <- total+num
      
      
    }
    
    # if(is.null(result[[1]])==0){
    #   temp <- length(result[[1]])
    #   snpid.result[total+(1:temp)] <- result[[1]]
    #   snpvalue.result[,total+(1:temp)] <- result[[2]]
    #   total <- temp+total
    # }
  }
  snpid.result <- snpid.result[1:total]
  snpvalue.result <- snpvalue.result[,1:total]
  
  
  conditional.snp.list.icog <- list(snpid.result,snpvalue.result)
  save(conditional.snp.list.icog,file=paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/conditional.snp.list.icog",i1,".Rdata"))
  
}


# ########## five snps only exists in ONCO array
# 
# extract.list.shared <- extract.list[1:total,]
# 
# save(extract.list.shared,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_list_shared.Rdata")
# extract.list.onco.only <- extract.list[(total+1):nrow(extract.list),]
# save(extract.list.onco.only,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_list_onco_only.Rdata")
# 
# 
# ##########
# idx.match <- match(extract.list.shared$SNP.ICOGS,snpid.result)
# snpid.result <- snpid.result[idx.match]
# all.equal(snpid.result,extract.list.shared$SNP.ICOGS)
# snpvalue.result <- snpvalue.result[,idx.match]
# extract.result <- list(snpid.result,snpvalue.result)
# colnames(snpvalue.result) <- snpid.result
# 
# 
# 
# extract.result <- list(snpid.result,snpvalue.result)
# save(extract.result,file="/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/ERPRHER2GRADE_casecase/result/extract_result_shared.Rdata")
# 
# 



















