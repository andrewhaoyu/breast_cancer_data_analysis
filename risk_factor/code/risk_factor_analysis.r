library(data.table)
library(bc2)
data <- fread("./data/dataset_montse_20180522.txt")
data.com <- fread("./data/concept_542-543-change-claude_bcac_pheno_v10_SPT_100217.txt")

data.com.new <- data.com[,c(2,40,48,56,10)]


new.


getZstandard()