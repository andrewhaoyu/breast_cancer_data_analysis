fine_mapping <- read.csv("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/data/fine_mapping_regions.csv",header=F,stringsAsFactors = F)

origin <- fine_mapping[,12]
origin <- gsub(".*:","",origin)
origin.split <- strsplit(origin,"_")
n <- nrow(fine_mapping)

start <- rep(0,n)
end <- rep(0,n)
for(i in 1:n){
  start[i] <- as.numeric(origin.split[[i]][1])
  end[i] <- as.numeric(origin.split[[i]][2])
}

fine_mapping <- cbind(fine_mapping,start,end)

write.csv(fine_mapping,file=paste0("/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/data/fine_mapping_regions.csv"),row.names =F ,col.names = F)
