results <- matrix(0,178,6)

for(i1 in 1:178){
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/known_SNPs/cluster_analysis/result",i1,".Rdata")) 
  results[i1,] <- dis.result
}

results.temp <- cbind(results[,1],results[,3])
cl <- kmeans(results,4)
plot(results[,1],results[,4],col=cl$cluster)





results.clean <- results
for(i in 1:6){
  idx <- which(results[,i]>=1.96^2)
  results.clean[idx,i] <- 1
  idx <- which(results[,i]<=1.96^2)
  results.clean[idx,i] <- 0
}





pca <- prcomp(results.clean)
plot(pca)
scores <- pca$x
plot(scores[,1],scores[,2])
pca.keep <- scores[,1:5]
cl <- kmeans(pca.keep,3)
plot(pca.keep[,1],pca.keep[,2],,col=cl$cluster)
cluster.vector <- cl$cluster
idx.1 <- which(cluster.vector==1)
