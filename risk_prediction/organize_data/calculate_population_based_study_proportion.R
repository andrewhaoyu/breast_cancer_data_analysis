#goal:calculate population based studies proportion
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
study_design_table = read.csv("./data/BCAC_study_design_simple.csv",header=T)

pop_study_cate = c("Population-based case-control study",
                   "Prospective cohort study: nested case-control",
                   "Population-based case-control",
                   "Nested case-control study",
                   "Case-control study, nested in a prospective cohort study",
                   "Population-based study of women <50 years",
                   "Cohort  study",
                   "Population-based",
                   "Population-based cohort study",
                   "Prospective cohort study: nested case-control study",
                   "Population-based prospective cohort study",
                   "Prospective Cohort Study (2003-2006) of women ages 35+ receiving screening mammography at Mayo Clinic and living in MN, IA, WI; nested case-control",
                   "Prospective cohort",
                   "Population-based case-control study, cohort study",
                   "Nested case control study from a population-based cohort",
                   "Nested case-control study of incident and prevalent cases within prospective cohort study",
                   "Population -based case-control study using prevalent cases with controls matched to cases on year of birth")
idx <- which(study_design_table$study_design%in%pop_study_cate)
pop_studies = study_design_table[idx,1]
write.table(pop_studies,file = "./data/pop_studies_defination.txt",quote = F, row.names = F, col.names = F)
data1_known = as.data.frame(fread("/data/zhangh24/breast_cancer_data_analysis/data/iCOGS_euro_v10_10232017.csv",header=T))

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
