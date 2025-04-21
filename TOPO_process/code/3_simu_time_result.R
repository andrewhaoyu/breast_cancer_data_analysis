#evaluate the relative time performance between TOPO and MTOP
setwd("/data/NC_BW/HZ_SF/simu_mtop/time")
result1 = fread("results_four_test100.csv")

mean(as.numeric(result1[3,2:4]/
       result1[1,2:4]))

result2 = fread("results_six_test100.csv")

mean(as.numeric(result2[3,2:4]/
                  result2[1,2:4]))

