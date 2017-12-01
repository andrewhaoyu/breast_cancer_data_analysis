data <- read.csv("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/data/fine_mapping_annotated_clean.csv",stringsAsFactors = F)
dim(data)
data <- data[1:178,]
head(data)
new.names <- data[,2]
old.names <- data[,1]
data.clean <- data[,2:ncol(data)]
idx.match <- match(old.names,new.names)
data.clean <- data.clean[idx.match,]
all.equal(data.clean[,1],old.names)
write.csv(data.clean,"/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/data/fine_mapping_annotated_clean.csv",quote=F,row.names =F)
