setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/discovery_SNP/enrichment_analysis/result/rediscoveryprojectnextsteps')
library(ggplot2)
###########220 celltypes results
celltypes220_lua <- read.table("220celltype.1.lumA.txt",
                               header=T)
celltypes220_lua$log10p <- -log10(celltypes220_lua$Enrichment_p)

celltypes220_TN <- read.table("220celltype.2.TN.txt",header=T)

celltypes220_TN$log10p <- -log10(celltypes220_TN$Enrichment_p)

# subtypes <- c(rep("Luminal A like",nrow(celltypes220_lua)),
#               rep("Triple negative",nrow(celltypes220_lua)))

################# create a data frame for gg plot2
n <- nrow(celltypes220_lua)
subtypes <- rep(c("Luminal A like","Triple negative"),n)
cell_types <- rep("c",2*n)
log10p <- rep(0,2*n)
for(i in 1:n){
  cell_types[2*i-1] <- paste0(as.character(celltypes220_lua$cell_type[i])," (",as.character(celltypes220_lua$mark[i]),")")
  cell_types[2*i] <- paste0(as.character(celltypes220_TN$cell_type[i])," (",as.character(celltypes220_TN$mark[i]), ")")
  log10p[2*i-1] <- celltypes220_lua$log10p[i]
  log10p[2*i] <-  celltypes220_TN$log10p[i]
}

celltype220_data <- data.frame(cell_types,log10p,subtypes)
colnames(celltype220_data) <- c("cell_types",
                                "log10p",
                                "subtypes")
celltype220_data1 <- celltype220_data[1:110,]
png("celltype220_data1.png", height = 32, width  = 28, res = 300)
ggplot(data=celltype220_data1,aes(x=cell_types,y=log10p))+
  geom_bar(stat = "identity",fill="steelblue")+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("-log10(p)") + 
  ggtitle("Enrichment analysis of 220 celltypes")  + 
  theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  
  geom_hline (yintercept = -log10(0.05/220), color = "red")+
  xlab("Cell types")+
  coord_flip()+
  facet_grid(.~subtypes)+
  theme(strip.text = element_text(face = "bold"))
dev.off()

celltype220_data2 <- celltype220_data[111:220,]
png("celltype220_data2.png", height = 32, width  = 28, res = 300)
ggplot(data=celltype220_data2,aes(x=cell_types,y=log10p))+
  geom_bar(stat = "identity",fill="steelblue")+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("-log10(p)") + 
  ggtitle("Enrichment analysis of 220 celltypes")  + 
  theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  
  geom_hline (yintercept = -log10(0.05/220), color = "red")+
  xlab("Cell types")+
  coord_flip()+
  facet_grid(.~subtypes)+
  theme(strip.text = element_text(face = "bold"))
dev.off()


celltype220_data3 <- celltype220_data[221:330,]
png("celltype220_data3.png", height = 32, width  = 28, res = 300)
ggplot(data=celltype220_data3,aes(x=cell_types,y=log10p))+
  geom_bar(stat = "identity",fill="steelblue")+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("-log10(p)") + 
  ggtitle("Enrichment analysis of 220 celltypes")  + 
  theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  
  geom_hline (yintercept = -log10(0.05/220), color = "red")+
  xlab("Cell types")+
  coord_flip()+
  facet_grid(.~subtypes)+
  theme(strip.text = element_text(face = "bold"))
dev.off()

celltype220_data4 <- celltype220_data[331:440,]
png("celltype220_data4.png", height = 32, width  = 28, res = 300)
ggplot(data=celltype220_data4,aes(x=cell_types,y=log10p))+
  geom_bar(stat = "identity",fill="steelblue")+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("-log10(p)") + 
  ggtitle("Enrichment analysis of 220 celltypes")  + 
  theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  
  geom_hline (yintercept = -log10(0.05/220), color = "red")+
  xlab("Cell types")+
  coord_flip()+
  facet_grid(.~subtypes)+
  theme(strip.text = element_text(face = "bold"))
dev.off()


#############baseline analysis
############################## create a data frame for gg plot2
baseline1.lua <- read.table("baseline.1.lumA.results.txt",header=T)
baseline2.TN <- read.table("baseline.2.TN.results.txt",header=T)

n <- nrow(baseline1.lua)
subtypes <- rep(c("Luminal A like","Triple negative"),n)
cell_types <- rep("c",2*n)
log10p <- rep(0,2*n)
for(i in 1:n){
  cell_types[2*i-1] <- paste0(as.character(baseline1.lua$Category[i]))
  cell_types[2*i] <- paste0(as.character(baseline2.TN$Category[i]))
  log10p[2*i-1] <- -log10(baseline1.lua$Enrichment_p[i])
  log10p[2*i] <-  -log10(baseline2.TN$Enrichment_p[i])
}

baseline_data <- data.frame(cell_types,log10p,subtypes)
colnames(baseline_data) <- c("Category",
                                "log10p",
                                "subtypes")
png("baseline_data.png", height = 32, width  = 28, res = 300)
ggplot(data=baseline_data,aes(x=Category,y=log10p))+
  geom_bar(stat = "identity",fill="steelblue")+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("-log10(p)") + 
  ggtitle("Enrichment analysis of ")  + 
  theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  
  geom_hline (yintercept = -log10(0.05/n), color = "red")+
  xlab("Cell types")+
  coord_flip()+
  facet_grid(.~subtypes)+
  theme(strip.text = element_text(face = "bold"))
dev.off()






