#test the batch effect for the 32 detected novel loci
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))
discovery_snp_paper <- read.csv("./data/discovery_snp_paper_order.csv")
head(discovery_snp_paper)
#test the difference of odds ratio between icogs and onco for first 22 snps from standard analysis
standard_snp <- discovery_snp_paper[1:22,]
standard_snp = standard_snp %>% 
  mutate(chr.pos = paste0(CHR,":",position))
data = data %>% 
  mutate(chr.pos=paste0(chr.Onco,":",Position.Onco))
standard_snp_all = left_join(standard_snp,
                             data,
                             by="chr.pos")
standard_snp_odds = standard_snp_all %>% 
  mutate(var.icogs = SE.iCOGs^2,
         var.onco = SE.Onco^2) %>% 
  select(SNP,CHR,position,beta.iCOGs,var.icogs,
         beta.Onco,var.onco)
TestDiff <- function(beta1,var1,beta2,var2){
   z= (beta1-beta2)/sqrt(var1+var2)
   p = 2*pnorm(-abs(z))
   return(p)
}
standard_snp_odds = standard_snp_odds %>% 
  mutate(diffP = TestDiff(beta.iCOGs,var.icogs,
                          beta.Onco,var.onco))
p1 = standard_snp_odds$diffP

#test the difference of odds ratio between icogs and onco for 8 snps from subtypes analysis
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ICOG/Intrinsic_subtypes/result/icog_result_shared_1p_082119.Rdata")
load("/spin1/users/zhangh24/breast_cancer_data_analysis/whole_genome_age/ONCO/intrinsic_subtypes/result/onco_result_shared_1p_082119.Rdata")

icog_logodds = icog_result_shared_1p %>%
  mutate(chr.pos = paste0(CHR,":",position))
  select(rs_id,chr.pos,11:40) %>% 
onco_logodds = onco_result_shared_1p %>% 
    mutate(chr.pos = paste0(CHR,":",position)) %>%
    select(rs_id,chr.pos,11:40)
  subtypes_snp <- discovery_snp_paper[23:32,]
  subtypes_snp = subtypes_snp %>% 
    mutate(chr.pos = paste0(CHR,":",position))
  subtypes_snp_all = left_join(subtypes_snp,
                               icog_logodds,
                               by="chr.pos")
  colnames(subtypes_snp_all)[15:19] <- paste0("icogs_beta",1:5)
  colnames(subtypes_snp_all)[20:44] <- paste0("icogs_var",1:25)
  subtypes_snp_all = left_join(subtypes_snp_all,
                               onco_logodds, 
                               by="chr.pos")
  colnames(subtypes_snp_all)[c(2:3)] <- c("CHR",
                                          "position")
  colnames(subtypes_snp_all)[c(50:54)] <- paste0("onco_beta",1:5)
  colnames(subtypes_snp_all)[c(55:79)] <- paste0("onco_var",1:25)
  subtypes_snp_odds <- subtypes_snp_all %>% 
    select(SNP,CHR,position,
           paste0("icogs_beta",1:5),
           paste0("icogs_var",1:25),
           paste0("onco_beta",1:5),
           paste0("onco_var",1:25))
  TestDiffvec <- function(beta1,var1,beta2,var2){
    betadiff= (beta1-beta2)
    vardiff = var1+var2
    chi2 = betadiff%*%solve(vardiff)%*%betadiff
    p = pchisq(chi2,df = length(betadiff))
    return(p)
  }
  n <- nrow(subtypes_snp_odds)
  p <- 5
p2 = rep(0,n)
for(i in 1:n){
  beta1 = as.numeric(subtypes_snp_odds[i,4:8])
  var1 = matrix(as.numeric(subtypes_snp_odds[i,9:33]),p,p)
  beta2 = as.numeric(subtypes_snp_odds[i,34:38])
  var2 = matrix(as.numeric(subtypes_snp_odds[i,39:63]),p,p)
  p2[i] = TestDiffvec(beta1,var1,beta2,var2)
}
