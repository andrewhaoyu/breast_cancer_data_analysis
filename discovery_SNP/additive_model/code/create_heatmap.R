#Goal: Create heatmap for discovery snp
setwd("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/")
data <- read.csv("./discovery_SNP/additive_model/result/additive_model_result.csv",header=T)
#install.packages("ComplexHeatmap")
library(tidyverse)
library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(colorspace)
library(circlize)
#library(xlsx)
library(readxl)
library(data.table)

GetZfun <- function(p){
  return(abs(qnorm(p/2)))
}
n.tu <- 4
snp <- data[,1]
new.result <- matrix(0,nrow(data),n.tu)
heter.p <- data[,22]
for(i in 1:n.tu){
  p <- data[,5+2*(i-1)]
  new.result[,i] <- GetZfun(p)
}
colnames(new.result) <- c("ER_z","PR_z",
                          "HER2_z","grade_z")

rownames(new.result) <- snp
png(filename = paste0("./discovery_SNP/additive_model/result/discovery_heatmap.png"),height = 20,width = 15,units ="cm",res = 300)
print(Heatmap(new.result,
              cluster_columns= FALSE,
              show_row_names = TRUE,
              row_names_side = "right",
              row_dend_side = "left",
              row_names_gp=gpar(cex=fontsize),
              row_dend_width = unit(3, "cm"),
              colorRamp2(c(1.281552, 1.959964, 4.891638), c("mistyrose3", "thistle1", "darkred")),
              heatmap_legend_param =list(color_bar =  "continuous", labels_gp = gpar(fontsize = 8), at = c(GetZfun(10^-6), GetZfun(10^-4), GetZfun(10^-3),GetZfun(0.05), GetZfun(0.2)), 
labels = c("10-6", "10-4", "10-3","0.05", "0.20"))
)
)
dev.off()

# ,
#         row_names_gp=gpar(cex=fontsize),
#         row_dend_width = unit(3, "cm"),
#         colorRamp2(c(1.281552, 1.959964, 6.10941), c("mistyrose3", "thistle1", "darkred")),
#         heatmap_legend_param =list(color_bar =  "continuous", labels_gp = gpar(fontsize = 8), at = c( 4.417173, 3.290527, 1.644854, 1.281552), 
#                                    labels = c("5.0x10-6", "5.0x10-4", "0.05", "0.20")))
# 
# heatmap.2(new.result)
# 
# new.result <- data.frame(snp,new.result,
#                          heter.p)
