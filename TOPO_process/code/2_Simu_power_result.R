#get the percentage of power
#increase for topo compared to alternative methods
#
setwd("/data/NC_BW/HZ_SF/simu_mtop/power")
library(data.table)
library(dplyr)
#load result four tumor characteristics
#column V1-V3 no heterogeneity setting (sample size:25,000; 50,000, and 100,000)
#column V4-V6 additive effects setting (sample size:25,000; 50,000, and 100,000)
#column V7-V9 interaction effects setting (sample size:25,000; 50,000, and 100,000)
result1 <- fread("./results2000_four.csv")
#load result six tumor characteristics
result2 <- fread("/data/BB_Bioinformatics/HZ/topo_result/global_power/results1980_six.csv")

result <- rbind(result1, result2)

#get the relative performance between MTOP and TOPO
method_name = result[,1]
result = result[,-1]
#connect four and six markers
as.numeric(c((result[6,4:6]-result[1,4:6]),
  (result[12,4:6]-result[7,4:6])))/
  as.numeric(c(result[1,4:6],
    result[7,4:6]))
#differences between topo and mtop under additive
mean(as.numeric(c((result[6,4:6]-result[1,4:6]),
                  (result[12,4:6]-result[7,4:6])))/
       as.numeric(c(result[1,4:6],
                    result[7,4:6]))
)

#differences between topo and standard logistic regression under additive
mean(as.numeric(c((result[6,4:6]-result[2,4:6]),
                  (result[12,4:6]-result[8,4:6])))/
       as.numeric(c(result[2,4:6],
                    result[8,4:6]))
)


#differences between topo and mtop under interaction
mean(as.numeric(c((result[6,7:9]-result[1,7:9]),
                  (result[12,7:9]-result[7,7:9])))/
       as.numeric(c(result[1,7:9],
                    result[7,7:9]))
)

#differences between topo and standard logistic regression under additive
x = (as.numeric(c((result[6,7:9]-result[2,7:9]),
                  (result[12,7:9]-result[8,7:9])))/
       as.numeric(c(result[2,7:9],
                    result[8,7:9])))
max_finite <- max(x[is.finite(x)])

# Replace Inf values with the maximum finite value
x[!is.finite(x)] <- max_finite

# Calculate the mean of the modified vector
mean_x <- mean(x)
