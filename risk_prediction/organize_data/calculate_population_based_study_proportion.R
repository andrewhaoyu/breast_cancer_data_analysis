#goal:calculate population based studies proportion
library(data.table)

data1_known = as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/iCOGS_euro_v10_10232017.csv",header=T))
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data2_known <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
library(dplyr)
data1_select = data1_known %>% 
  select(SG_ID,study,Behaviour1,age) %>% 
  rename(ID = SG_ID)
data2_select = data2_known %>% 
  select(Onc_ID,study,Behaviour1,age)%>% 
  rename(ID = Onc_ID)

data_com = rbind(data1_select, data2_select) %>% 
  filter(age!=888)

pop_studies = c("ABCFS", "BCEES", "BCINIS", "BREOGAN", "CAMA", "CBCS","CECILE","CGPS","ESTHER","GENICA",
                "KOHBRA","LAABC","LIFEPOOL","MASTOS","MISS","MTLGEBCS","NBHS","PBCS","PROCAS",
                "SASBAC","SBCGS","SEARCH","UCIBCS","US3SS","USRT")
n_case = length(which(data_com$Behaviour1==1))
n_control = length(which(data_com$Behaviour1==0))
data_pop = data_com %>% 
  filter(study%in%pop_studies)
#GESBC Population-based study of women <50 years
#HCSC Population-based prospective clinical cohort
#KARBAC Population and hospital-based cases; geographically matched controls
#KBCP Population-based prospective clinical cohort
#NC-BCFR Population-based recruitment of families; family-based cohort; population-based controls for subset of cases
#OFBCR Population-based familial case-control study 
#SMC Nested case control study from a population-based cohort
n_case_pop =  length(which(data_pop$Behaviour1==1))
n_control_pop =  length(which(data_pop$Behaviour1==0))
pop_pro = nrow(data_pop)/nrow(data_com)
case_pro = n_case_pop/n_case
control_pro = n_control_pop/n_control

result = data.frame(overall = pop_pro, case = case_pro, control_pro)
write.csv(result, file = "/data/zhangh24/breast_cancer_data_analysis/data/population_based_study_pro.csv",
          quote = F, row.names = F)

age = data_com$age
age_control = data_com %>% filter(Behaviour1==0) %>% select(age) 
age_control = age_control$age
age_case = data_com %>% filter(Behaviour1==1) %>% select(age)
age_case = age_case$age
age_mean = mean(age)
age_sd = sd(age)
age_mean_control = mean(age_control)
age_sd_control = sd(age_control)
age_mean_case = mean(age_case)
age_sd_case = sd(age_case)
age_med = median(age)
age_iqr = IQR(age)
age_med_control = median(age_control)
age_iqr_control = IQR(age_control)
age_med_case = median(age_case)
age_iqr_case = IQR(age_case)

result = data.frame(mean = c(age_mean,age_mean_case,age_mean_control),
                    sd = c(age_sd,age_sd_case,age_sd_control),
                    median = c(age_med,age_med_case,age_med_control),
                    iqr = c(age_iqr,age_iqr_case,age_iqr_control))
row.names(result) = c("overall", "case", "control")
write.csv(result, file = "/data/zhangh24/breast_cancer_data_analysis/data/age_distribution.csv",
          quote = F, row.names = T)
