setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/discovery_SNP/enrichment_analysis/result/rediscoveryprojectnextsteps')
library(ggplot2)
library(reshape)
###########220 celltypes results

celltypes220_lua <- read.table("220celltype.1.lumA.txt",
                               header=T)
celltypes220_lua$log10p <- -log10(celltypes220_lua$Enrichment_p)

celltypes220_TN <- read.table("220celltype.2.TN.txt",header=T)

celltypes220_TN$log10p <- -log10(celltypes220_TN$Enrichment_p)

# subtypes <- c(rep("Luminal A like",nrow(celltypes220_lua)),
#               rep("Triple negative",nrow(celltypes220_lua)))
#################create a plot for enrichment score
n <- nrow(celltypes220_lua)
cell_type <- paste0(as.character(celltypes220_lua$cell_type), "_",
                    as.character(celltypes220_lua$mark))
cell_type[137] <- "Fetal_brain_(H3K4me3)_2"
cell_type[60] <- "Pancreatic_islets_H3K4me1_2"
cell_type[138] <- "Pancreatic_islets_H3K4me3_2"
luminal_enrichment_score <- celltypes220_lua$Coefficient_z.score
triple_neg_enrichment_score <- celltypes220_TN$Coefficient_z.score
celltype220_data <- data.frame(luminal_enrichment_score,
                               triple_neg_enrichment_score)
rownames(celltype220_data) <- cell_type
colnames(celltype220_data) <- c("Luminal A like",
                                "Triple negative")
celltype220_data <- as.matrix(celltype220_data)
#plot(celltype220_data[,1],celltype220_data[,2])
celltype220_data1 <- celltype220_data[c(55,56,57,58,133,134,135,181),]
library(pheatmap)
library(RColorBrewer)
library(viridis)
paletteLength <- 20
my.color <- colorRampPalette(c("dodgerblue4", "white", "#c0392b"))(paletteLength)
myBreaks <- c(seq(min(celltype220_data1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(celltype220_data1)/paletteLength, max(celltype220_data1), length.out=floor(paletteLength/2)))
png("brest_cell_z.png", height = 32, width  = 28,units="cm" ,res = 300)
pheatmap(
  mat               = celltype220_data1,
  color             = my.color,
  breaks           = myBreaks,
  cluster_rows=FALSE,
  cluster_cols=FALSE,
  # border_color      = NA,
  # show_colnames     = FALSE,
  # show_rownames     = FALSE,
  # annotation_col    = mat_col,
  # annotation_colors = mat_colors,
  # drop_levels       = TRUE,
  fontsize          = 14
  # main              = "Default Heatmap"
)
dev.off()

###############generate the plot for H3K27ac
idx <- which(celltypes220_lua$mark=="H3K27ac")

celltype220_data2 <- celltype220_data[idx,]
# library(pheatmap)
# library(RColorBrewer)
# library(viridis)
paletteLength <- 50
my.color <- colorRampPalette(c("dodgerblue4", "white", "#c0392b"))(paletteLength)
myBreaks <- c(seq(min(celltype220_data2), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(celltype220_data2)/paletteLength, max(celltype220_data1), length.out=floor(paletteLength/2)))
png("H3K27_z.png", height = 32, width  = 28, units="cm",res = 300)
pheatmap(
  mat               = celltype220_data2,
  color             = my.color,
  breaks           = myBreaks,
  cluster_rows=FALSE,
  cluster_cols=FALSE,
  # border_color      = NA,
  # show_colnames     = FALSE,
  # show_rownames     = FALSE,
  # annotation_col    = mat_col,
  # annotation_colors = mat_colors,
  # drop_levels       = TRUE,
  fontsize          = 12
  # main              = "Default Heatmap"
)
dev.off()






celltype220_data1.m <- melt(celltype220_data1)
base_size <- 9
ggplot(celltype220_data1.m,aes(X2,X1))+
  geom_tile(aes(fill=value),colour="white")+
  scale_fill_gradient(low="white",
                      high="steelblue")+
  theme_grey(base_size = base_size) + 
  labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) 
  # element(legend.position = "none", 
  #      axis.ticks = theme_blank(), 
  #      axis.text.x = theme_text(size = base_size *0.8, angle = 330, hjust = 0, colour = "grey50"))

#heatmap(celltype220_data1,scale = "column")




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






