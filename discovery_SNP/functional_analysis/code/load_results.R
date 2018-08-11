setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/discovery_SNP/functional_analysis/result/rediscoveryprojectnextsteps')
celltypes220_lua <- read.table("220celltypes.1.lumA.txt",
                               header=T)
celltypes220_TN <- read.table("220celltypes.2.TN.txt",header=T)

baseline1.lua <- read.table("baseline.1.lumA.results",header=T)

baseline2.TN <- read.table("baseline.2.TN.results",header=T)

eQTL.lua <- read.table("eQTL.1.lumA.txt",header=T)

eQTL.TN <-  read.table("eQTL.2.TN.txt",header=T)

