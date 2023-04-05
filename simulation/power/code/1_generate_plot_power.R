data = read.csv("./simulation/power/result/TOPO_results/power_four_makrer.csv",header = F)
data[1,1] = "sample_size"
colnames(data)[2:10] = c(rep("No Hetergeneity",3),
                          rep("Additive Effect",3),
                          rep("Additive + Interaction Effects",3))
method = rep("c", 100)
sample_size = rep(0, 100)
ga_setting = rep("c", 100)
power = rep(0, 100)
temp = 1
for(i in 2:nrow(data)){
  for(k in 2:ncol(data)){
    sample_size[temp] = data[1, k]
    power[temp] = data[i,k]
    method[temp] = data[i,1]
    ga_setting[temp] = colnames(data)[k]
    temp = temp + 1
  }
}
method = method[1:(temp-1)]
sample_size = sample_size[1:(temp-1)]
ga_setting = ga_setting[1:(temp-1)]
power = power[1:(temp-1)]
source("/Users/zhangh24/GoogleDrive/multi_ethnic//code/LD_simulation_large/theme_Publication.R")
library(tidyverse)
plot_data = data.frame(method, sample_size, ga_setting, power)

plot_data = plot_data %>% mutate(
  method_new = case_when(
    method%in% "poly" ~ "Polytomous model",
    method%in% "standard" ~ "Standard logistic regression",
    method%in% "topo" ~ "TOPO",
    method%in% "MTOP" ~ "MTOP"
  ),
  sample_size = as.character(sample_size),
  sample_size_new = case_when(
    sample_size %in% "25000" ~ "25000",
    sample_size %in% "50000" ~ "50000",
    sample_size %in% "1e+05" ~ "100000"
  ),
  sample_size_new = factor(sample_size_new,
                           levels = c("25000",
                                      "50000",
                                      "100000")),
  method_new = factor(method_new, 
                      levels = c("Standard logistic regression",
                                 "MTOP",
                                 "TOPO",
                                 "Polytomous model")),
  ga_setting = factor(ga_setting)
) %>% 
select(method_new, sample_size_new, ga_setting, power) %>% 
  rename(method = method_new,
         sample_size = sample_size_new)
plot_data_sub = plot_data %>% 
  filter(ga_setting %in% "No Hetergeneity")
p <- ggplot(plot_data_sub)+
  geom_bar(aes(x = sample_size,y = power,fill=method),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab(expression(bold(paste("Power"))))+
  theme_Publication()
png(filename = "./simulation/power/result/TOPO_results/no_hetergeneity.png",
    width = 12, height = 8, unit = "in", res = 300)
print(p)
dev.off()
plot_data_sub = plot_data %>% 
  filter(ga_setting %in% "Additive Effect")
p <- ggplot(plot_data_sub)+
  geom_bar(aes(x = sample_size,y = power,fill=method),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab(expression(bold(paste("Power"))))+
  theme_Publication()
png(filename = "./simulation/power/result/TOPO_results/additive_effect.png",
    width = 12, height = 8, unit = "in", res = 300)
print(p)
dev.off()
plot_data_sub = plot_data %>% 
  filter(ga_setting %in% "Additive + Interaction Effects")

p <- ggplot(plot_data_sub)+
  geom_bar(aes(x = sample_size,y = power,fill=method),
           position = position_dodge(),
           stat = "identity")+
  theme_Publication()+
  ylab(expression(bold(paste("Power"))))+
  theme_Publication()
png(filename = "./simulation/power/result/TOPO_results/interaction_effect.png",
    width = 12, height = 8, unit = "in", res = 300)
print(p)
dev.off()
