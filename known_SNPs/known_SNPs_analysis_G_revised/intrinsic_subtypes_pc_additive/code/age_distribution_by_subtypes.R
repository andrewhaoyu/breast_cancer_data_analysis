library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")
data1 <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data2 <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
library(dplyr)



data1_select = data1 %>% 
  select(BCAC_ID,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1,age) %>% 
  rename(case = Behaviour1,
         ER = ER_status1,
         PR = PR_status1,
         HER2 = HER2_status1,
         Grade = Grade1)
data2_select = data2 %>% 
  select(BCAC_ID,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1,age) %>% 
  rename(case = Behaviour1,
         ER = ER_status1,
         PR = PR_status1,
         HER2 = HER2_status1,
         Grade = Grade1)
data_select = rbind(data1_select,data2_select) 
data_select_filter = data_select %>% filter(age!=888)
data_select_filter %>% group_by(case) %>% 
  summarize(mean(age))
mean(data_select_filter$age)
# %>%
#   filter(age!=888)

#control
control = data_select %>% filter(case==0)
quantile(control$age)
#luminal A
idx.1 <- which((data_select$ER==1|data_select$PR==1)
               &data_select$HER2==0
               &(data_select$Grade==1|data_select$Grade==2))
lua = data_select[idx.1,]



idx.1 <- which((data1_select$ER==1|data1_select$PR==1)
               &data1_select$HER2==0
               &(data1_select$Grade==1|data1_select$Grade==2))
lua = data1_select[idx.1,]

idx.1 <- which((data2_select$ER==1|data2_select$PR==1)
               &data2_select$HER2==0
               &(data2_select$Grade==1|data2_select$Grade==2))
lua = data2_select[idx.1,]
dim(lua)


#quantile(lua$age)
table(cut(lua$age,breaks = c(0,40,50,60,100,1000),right = F))
#luminal B
idx.2 <- which((data_select$ER==1|data_select$PR==1)
               &data_select$HER2==1)
lub = data_select[idx.2,]
table(cut(lub$age,breaks = c(0,40,50,60,100,1000),right = F))
#HR+_HER2-_highgrade
idx.3 <- which((data_select$ER==1|data_select$PR==1)
               &data_select$HER2==0
               &data_select$Grade==3)
lub_hern = data_select[idx.3,]
table(cut(lub_hern$age,breaks = c(0,40,50,60,100,1000),right = F))
#for third subtype HR-_HER2+
idx.4 <- which(data_select$ER==0&data_select$PR==0
               &data_select$HER2==1)
hern = data_select[idx.4,]
table(cut(hern$age,breaks = c(0,40,50,60,100,1000),right = F))
#for third subtype HR-_HER2-
idx.5 <- which(data_select$ER==0&data_select$PR==0
               &data_select$HER2==0)
TN = data_select[idx.5,]
table(cut(TN$age,breaks = c(0,40,50,60,100,1000),right = F))

quantile_table = rbind(table(cut(control$age,breaks = c(0,40,50,60,100,1000),right = F)),
                       table(cut(lua$age,breaks = c(0,40,50,60,100,1000),right = F)),
                       table(cut(lub$age,breaks = c(0,40,50,60,100,1000),right = F)),
                       table(cut(lub_hern$age,breaks = c(0,40,50,60,100,1000),right = F)),
                       table(cut(hern$age,breaks = c(0,40,50,60,100,1000),right = F)),
                       table(cut(TN$age,breaks = c(0,40,50,60,100,1000),right = F)))
row.names(quantile_table) = c("Control",
                              "HR+_HER2-_lowgrade",                              
                              "HR+_HER2+",
                              "HR+_HER2-_highgrade",
                              "HR-_HER2+",
                              "HR-_HER2-")
colnames(quantile_table) = c("<40","40-49","50-59","60+","unknown")
write.csv(quantile_table,file = "/data/zhangh24/breast_cancer_data_analysis/data/age_distribution.csv")
data1_select = data1 %>% 
  select(SG_ID,BCAC_ID,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1,age)
  write.csv(data1_select,file = "/data/zhangh24/breast_cancer_data_analysis/data/icogs_id.csv",row.names = F)
  data2_select = data2 %>% 
    select(Onc_ID,BCAC_ID,Behaviour1,ER_status1,PR_status1,HER2_status1,Grade1,age)
    write.csv(data2_select,file = "/data/zhangh24/breast_cancer_data_analysis/data/onco_id.csv",row.names = F)
 
