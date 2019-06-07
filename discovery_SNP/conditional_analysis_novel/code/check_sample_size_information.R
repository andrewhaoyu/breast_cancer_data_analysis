#Goal: Check sample size information
data1 <- as.data.frame(fread("./data/icogs_overall.csv",header=T))
#filter icogs greece and Belarus
library(dplyr)
filter1 = data1 %>% 
  filter(StudyCountry=="Belarus"|
           StudyCountry=="Greece") %>% 
  select(Behaviour1,StudyCountry)
table(filter1$Behaviour1)
A_filter1 = data1 %>% 
  filter(StudyCountry!="Belarus"&
           StudyCountry!="Greece"&
           Behaviour1!=2&
           Behaviour1!=888) %>% 
  select(Behaviour1,StudyCountry)

table(A_filter1$Behaviour1)

data2 <- fread("./data/oncoarray_overall.csv",header=T)
filter2 = data2 %>% 
  filter(StudyCountry=="Belarus") %>% 
  select(Behaviour1,StudyCountry)
table(filter2$Behaviour1)
A_filter2 = data2 %>% 
  filter(StudyCountry!="Belarus"&
           Behaviour1!=2&
           Behaviour1!=888)%>% 
  select(Behaviour1,StudyCountry)
table(A_filter2$Behaviour1)

data1.new <- fread("./data/iCOGS_euro_v10_10232017.csv",header=T)
data2.new <- fread("./data/Onco_euro_v10_10232017.csv",header=T)
table(data2.new$Behaviour1)
