setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis/discovery_SNP/enrichment_analysis/result/rediscoveryprojectnextsteps')
library(ggplot2)
library(reshape)

##########test difference

TestDiff <- function(x1,s1,x2,s2){
  z = (x1-x2)/sqrt(s1^2+s2^2)
  return(2*pnorm(abs(z),lower.tail = F))
}

##########baseline results main
baseline_result <- read.csv("baseline_results_main.csv")
png("baseline_main_results.png", height = 16, width  = 19.04, res = 300, units = "cm")
low.95 <- baseline_result$Enrichment - baseline_result$Enrichment_std_error
high.95 <- baseline_result$Enrichment + baseline_result$Enrichment_std_error
baseline_result$high.95 <- high.95
baseline_result$low.95 <- low.95
ggplot(data=baseline_result,aes(x=Annotation,y=Enrichment))+
  geom_bar(stat = "identity",
           position = "dodge",
           aes(fill=Subtypes)
           )+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("Enrichment") + 
  #ggtitle("Enrichment analysis of 220 celltypes")  + 
  theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold",hjust = 0.5,angle = 90),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  geom_hline(yintercept = 1, linetype="dashed", color = "blue")+
  coord_flip()+
  scale_y_continuous(breaks=c(-5,0,1,5,10,15,20))+
  geom_errorbar(aes(ymax = high.95,
                      ymin = low.95,
                      shape=Subtypes),
                  position= "dodge")+
    scale_color_manual(values=c("black","black"))

  #geom_hline (yintercept = -log10(0.05/220), color = "red")+
  #xlab("Cell types")+
  
  #facet_grid(.~subtypes)+
#  theme(strip.text = element_text(face = "bold"))
dev.off()


baseline_result <- read.csv("baseline_results_main.csv")
png("baseline_main_results_95.png", height = 16, width  = 19.04, res = 300, units = "cm")
low.95 <- baseline_result$Enrichment - baseline_result$Enrichment_std_error
high.95 <- baseline_result$Enrichment + baseline_result$Enrichment_std_error
baseline_result$high.95 <- high.95
baseline_result$low.95 <- low.95
ggplot(data=baseline_result,aes(x=Annotation,y=Enrichment))+
  geom_bar(stat = "identity",
           position = "dodge",
           aes(fill=Subtypes)
  )+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("Enrichment") + 
  #ggtitle("Enrichment analysis of 220 celltypes")  + 
  theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold",hjust = 0.5,angle = 90),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  #geom_hline(yintercept = 1, linetype="dashed", color = "blue")+
  coord_flip()+
  scale_y_continuous(breaks=c(-5,0,5,10,15,20))+
 geom_errorbar(aes(ymax = high.95, 
                     ymin = low.95,
                     shape=Subtypes),
                 position= "dodge")+
   scale_color_manual(values=c("black","black"))
# 
#geom_hline (yintercept = -log10(0.05/220), color = "red")+
#xlab("Cell types")+

#facet_grid(.~subtypes)+
#  theme(strip.text = element_text(face = "bold"))
dev.off()





###########baseline results with 500kb extension
baseline_result <- read.csv("baseline_500kb.csv")

low.95 <- baseline_result$Enrichment - baseline_result$Enrichment_std_error
high.95 <- baseline_result$Enrichment + baseline_result$Enrichment_std_error
baseline_result$high.95 <- high.95
baseline_result$low.95 <- low.95

png("baseline_results_500bp.png",height = 16, width  = 19.04, res = 300, units = "cm" )
ggplot(data=baseline_result,aes(x=Annotation,y=Enrichment))+
  geom_bar(stat = "identity",
           position = "dodge",
           aes(fill=Subtypes)
  )+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("Enrichment") + 
  #ggtitle("Enrichment analysis of 220 celltypes")  + 
  theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold",hjust = 0.5,angle = 0),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
 geom_hline(yintercept = 1, linetype="dashed", color = "blue")+
  coord_flip()+
  scale_y_continuous(breaks=c(-5,0,1,5,10,15,20))+
geom_errorbar(aes(ymax = high.95,
                    ymin = low.95,
                    shape=Subtypes),
                position= "dodge")+
  scale_color_manual(values=c("black","black"))
dev.off()

png("baseline_results_500bp_95.png",height = 16, width  = 19.04, res = 300, units = "cm" )
ggplot(data=baseline_result,aes(x=Annotation,y=Enrichment))+
  geom_bar(stat = "identity",
           position = "dodge",
           aes(fill=Subtypes)
  )+
  theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("Enrichment") + 
  #ggtitle("Enrichment analysis of 220 celltypes")  + 
  theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold",hjust = 0.5,angle = 0),
        axis.text.y = element_text(face = "bold.italic"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))+
  #geom_hline(yintercept = 1, linetype="dashed", color = "blue")+
  coord_flip()+
  scale_y_continuous(breaks=c(-5,0,5,10,15,20))+
geom_errorbar(aes(ymax = high.95,
                     ymin = low.95,
                     shape=Subtypes),
                 position= "dodge")+
   scale_color_manual(values=c("black","black"))
dev.off()




#############baseline analysis
############################## create a data frame for gg plot2
baseline1.lua <- read.table("baseline.1.lumA.results",header=T)
baseline2.TN <- read.table("baseline.2.TN.results",header=T)

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
# png("baseline_data.png", height = 32, width  = 28, res = 300,
#     units = "cm")
# ggplot(data=baseline_data,aes(x=Category,y=log10p))+
#   geom_bar(stat = "identity",fill="steelblue")+
#   theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("-log10(p)") + 
#   ggtitle("Enrichment analysis of ")  + 
#   theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold"),
#         axis.text.y = element_text(face = "bold.italic"),
#         axis.title.x = element_text(face = "bold"),
#         axis.title.y = element_text(face = "bold"))+
#   
#   geom_hline (yintercept = -log10(0.05/n), color = "red")+
#   xlab("Cell types")+
#   coord_flip()+
#   facet_grid(.~subtypes)+
#   theme(strip.text = element_text(face = "bold"))
# dev.off()
# 


############analysze the baseline results

idx1 <- which(baseline1.lua$Enrichment_p<=0.05/52)
baseline1.lua[idx1,]
idx2 <- which(baseline2.TN$Enrichment_p<=0.05/52)
baseline2.TN[idx2,]
baseline2.TN[order(baseline2.TN$Enrichment_p),]
TestDiff(baseline1.lua$Enrichment,
         baseline1.lua$Enrichment_std_error,
         baseline2.TN$Enrichment,
         baseline2.TN$Enrichment_std_error)
p.value <- TestDiff(baseline1.lua$Coefficient,
                    baseline1.lua$Coefficient_std_error,
                    baseline2.TN$Coefficient,
                    baseline2.TN$Coefficient_std_error)
write.csv(p.value,"baseline_test_difference.csv")
idx <- which.min(p.value)
order(p.value)
baseline1.lua[idx,]
baseline1.lua[46,]
baseline1.lua$Enrichment-baseline2.TN$Enrichment
sqrt(baseline1.lua$Enrichment_std_error^2+
       baseline2.TN$Enrichment_std_error^2)
write.csv(baseline2.TN,file = "baseline_TN.csv")
write.csv(baseline1.lua,file = "baseline_Lua.csv")

















###########220 celltypes results
###########breast related tissues: 55, 56, 57, 58, 133, 134, 135, 181, 275, 276, 277, 278, 353, 354, 355, 401
###########
celltypes220 <- read.csv("celltypes_result.csv",header=T,sep ="\t")
celltypes_breast <- celltypes220[c(55, 56, 57, 58, 133, 134, 135, 181, 275, 276, 277, 278, 353, 354, 355, 401),]
celltypes_breast$Annotation <- paste0(celltypes_breast$cell_type,"_",
                                      celltypes_breast$mark)


# ggplot(data=celltypes_breast,aes(x=Annotation,y=Enrichment))+
#   geom_bar(stat = "identity",
#            position = "dodge",
#            aes(fill=Subtypes)
#   )+
#   theme_minimal()+  theme(axis.text.x = element_text(angle = 0, hjust = 1)) + ylab("Enrichment") + 
#   #ggtitle("Enrichment analysis of 220 celltypes")  + 
#   theme(text = element_text(size=10),plot.title = element_text(hjust = 0.5,face = "bold"),axis.text.x = element_text(face = "bold",hjust = 0.5,angle = 90),
#         axis.text.y = element_text(face = "bold.italic"),
#         axis.title.x = element_text(face = "bold"),
#         axis.title.y = element_text(face = "bold"))+
#   geom_hline(yintercept = 1, linetype="dashed", color = "blue")+
#   coord_flip()+
#   scale_y_continuous(breaks=c(-5,0,1,5,10,15,20))
# 
# 


celltypes220_lua <- read.csv("220celltype.1.lumA.txt",
                                header=T,sep= "\t")
# write.csv(celltypes220_lua,file = "celltype2220_lumA.csv")
 celltypes220_lua$log10p <- -log10(celltypes220_lua$Enrichment_p)
# 
 celltypes220_TN <- read.table("220celltype.2.TN.txt",header=T)
# write.csv(celltypes220_TN,file = "celltypes220_TN.csv")



idx1 <- which(celltypes220_lua$Enrichment_p<=0.05/220)
celltypes220_lua[idx1,]
idx2 <- which(celltypes220_TN$Enrichment_p<=0.05/220)
celltypes220_TN[idx2,]
library(xlsx)
write.csv(celltypes220_lua[idx1,],file = "sig_220_luA.csv")
result.celltypes <- rbind(celltypes220_lua[idx1,],
                         celltypes220_TN[idx2,])
##################test the difference between the 220 celltypes
p.value.220 <- TestDiff(celltypes220_lua$Coefficient,celltypes220_lua$Coefficient_std_error,celltypes220_TN$Coefficient,celltypes220_TN$Coefficient_std_error)
idx <- which(celltypes.diff<=0.05/220)
print(idx)



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
big_cell_type <- read.csv("10_celltypes_group.csv")
celltype220_data <- data.frame(luminal_enrichment_score,
                               triple_neg_enrichment_score
                               )
cell_type_10 <- data.frame(big_cell_type$Big_cell_type)
colnames(big_cell_type) <- "Cell type"
rownames(celltype220_data) <- cell_type
colnames(celltype220_data) <- c("HR+, HER2-, low grade",
                                "Triple negative"
                                )
celltype220_data <- as.matrix(celltype220_data)
#plot(celltype220_data[,1],celltype220_data[,2])
############plot for breast related tissues
celltype220_data1 <- celltype220_data[c(55,56,57,58,133,134,135,181),]
p.value.220[c(55,56,57,58,133,134,135,181)]
library(pheatmap)
library(RColorBrewer)
library(viridis)
paletteLength <- 20
my.color <- colorRampPalette(c("dodgerblue4", "white", "#c0392b"))(paletteLength)
myBreaks <- c(seq(min(celltype220_data1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(celltype220_data1)/paletteLength, max(celltype220_data1), length.out=floor(paletteLength/2)))
png("breast_cell_z.png", height = 16, width  = 18.5,units="cm" ,res = 300)
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


#############color setting
annotation_colors_full = list(group = brewer.pal(10, "Paired"))
names(annotation_colors_full$group) <- c("Adipose",
                                    "Other",
                                    "CNS",
                                    "Immune",
                                    "Gastrointestinal",
                                    "Cardiovascular",
                                    "Adrenal/Pancreas",
                                    "Skeletal muscle",
                                    "Breast",
                                    "Skin")













###############generate the plot for H3K27ac
idx <- which(celltypes220_lua$mark=="H3K27ac")

celltype220_data_H3K27ac <- celltype220_data[idx,]
row.names(celltype220_data_H3K27ac) <- as.character(celltypes220_lua$cell_type[idx])
group_H3K27ac <- data.frame(group = as.character(cell_type_10[idx,]),stringsAsFactors = F)
row.names(group_H3K27ac) <- row.names(celltype220_data_H3K27ac) 
group_H3K27ac <- data.frame(group_H3K27ac)
colnames(group_H3K27ac) <- "Cell type"
# library(pheatmap)
# library(RColorBrewer)
# library(viridis)
paletteLength <- 50
annotation_colors = annotation_colors_full
annotation_colors$group <- annotation_colors_full$group[names(annotation_colors_full$group)%in%as.character(unique(cell_type_10[idx,]))]
names(annotation_colors) <- "Cell type"
my.color <- colorRampPalette(c("dodgerblue4", "white", "#c0392b"))(paletteLength)
myBreaks <- c(seq(min(celltype220_data_H3K27ac), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(celltype220_data_H3K27ac)/paletteLength, max(celltype220_data_H3K27ac), length.out=floor(paletteLength/2)))
png("H3K27_z.png", height = 1290, width  = 1017)
#par(font.axis = 2)
pheatmap(
  mat               = celltype220_data_H3K27ac,
  color             = my.color,
  breaks           = myBreaks,
  cluster_rows=T,
  cluster_cols=F,
  annotation_row = group_H3K27ac,
  annotation_colors = annotation_colors,
  # border_color      = NA,
  # show_colnames     = FALSE,
  # show_rownames     = FALSE,
  # annotation_col    = mat_col,
  # annotation_colors = mat_colors,
  # drop_levels       = TRUE,
  fontsize          = 13,
  cellwidth=220, cellheight=16
  # main              = "Default Heatmap"
)
#axis(side = 2, font = 2)
dev.off()

###############generate the plot for H3K4me1
idx <- which(celltypes220_lua$mark=="H3K4me1")

celltype220_data_H3K4me1 <- celltype220_data[idx,]
row.names(celltype220_data_H3K4me1) <- as.character(big_cell_type[idx,1])
group_H3K4me1 <- data.frame(group = as.character(cell_type_10[idx,]),stringsAsFactors = F)
row.names(group_H3K4me1) <- row.names(celltype220_data_H3K4me1) 
group_H3K4me1 <- data.frame(group_H3K4me1)
colnames(group_H3K4me1) <- "Cell type"
# library(pheatmap)
# library(RColorBrewer)
# library(viridis)
paletteLength <- 50
annotation_colors = annotation_colors_full
annotation_colors$group <- annotation_colors_full$group[names(annotation_colors_full$group)%in%as.character(unique(cell_type_10[idx,]))]
names(annotation_colors) <- "Cell type"
my.color <- colorRampPalette(c("dodgerblue4", "white", "#c0392b"))(paletteLength)
myBreaks <- c(seq(min(celltype220_data_H3K4me1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(celltype220_data_H3K4me1)/paletteLength, max(celltype220_data_H3K4me1), length.out=floor(paletteLength/2)))
png("H3K4me1_z.png", height = 1290, width  = 1195)
pheatmap(
  mat               = celltype220_data_H3K4me1,
  color             = my.color,
  breaks           = myBreaks,
  cluster_rows=T,
  cluster_cols=F,
  annotation_row = group_H3K4me1,
  annotation_colors = annotation_colors,
  # border_color      = NA,
  # show_colnames     = FALSE,
  # show_rownames     = FALSE,
  # annotation_col    = mat_col,
  # annotation_colors = mat_colors,
  # drop_levels       = TRUE,
  fontsize          = 11,
  cellwidth=180, cellheight=10
  # main              = "Default Heatmap"
)
dev.off()





###############generate the plot for H3K4me3
idx <- which(celltypes220_lua$mark=="H3K4me3")

celltype220_data_H3K4me3 <- celltype220_data[idx,]
row.names(celltype220_data_H3K4me3) <- as.character(big_cell_type[idx,1])
group_H3K4me3 <- data.frame(group = as.character(cell_type_10[idx,]),stringsAsFactors = F)
row.names(group_H3K4me3) <- row.names(celltype220_data_H3K4me3) 
group_H3K4me3 <- data.frame(group_H3K4me3)
colnames(group_H3K4me3) <- "Cell type"
# library(pheatmap)
# library(RColorBrewer)
# library(viridis)
paletteLength <- 50
annotation_colors = annotation_colors_full
annotation_colors$group <- annotation_colors_full$group[names(annotation_colors_full$group)%in%as.character(unique(cell_type_10[idx,]))]
names(annotation_colors) <- "Cell type"
my.color <- colorRampPalette(c("dodgerblue4", "white", "#c0392b"))(paletteLength)
myBreaks <- c(seq(min(celltype220_data_H3K4me3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(celltype220_data_H3K4me3)/paletteLength, max(celltype220_data_H3K4me3), length.out=floor(paletteLength/2)))
png("H3K4me3_z.png", height = 1290, width  = 1195)
pheatmap(
  mat               = celltype220_data_H3K4me3,
  color             = my.color,
  breaks           = myBreaks,
  cluster_rows=T,
  cluster_cols=F,
  annotation_row = group_H3K4me3,
  annotation_colors = annotation_colors,
  # border_color      = NA,
  # show_colnames     = FALSE,
  # show_rownames     = FALSE,
  # annotation_col    = mat_col,
  # annotation_colors = mat_colors,
  # drop_levels       = TRUE,
  fontsize          = 11,
  cellwidth=180, cellheight=10
  # main              = "Default Heatmap"
)
dev.off()



###############generate the plot for H3K9ac
idx <- which(celltypes220_lua$mark=="H3K9ac")

celltype220_data_H3K9ac <- celltype220_data[idx,]
row.names(celltype220_data_H3K9ac) <- as.character(big_cell_type[idx,1])
group_H3K9ac <- data.frame(group = as.character(cell_type_10[idx,]),stringsAsFactors = F)
row.names(group_H3K9ac) <- row.names(celltype220_data_H3K9ac) 
group_H3K9ac <- data.frame(group_H3K9ac)
colnames(group_H3K9ac) <- "Cell type"
# library(pheatmap)
# library(RColorBrewer)
# library(viridis)
paletteLength <- 50
annotation_colors = annotation_colors_full
annotation_colors$group <- annotation_colors_full$group[names(annotation_colors_full$group)%in%as.character(unique(cell_type_10[idx,]))]
names(annotation_colors) <- "Cell type"
my.color <- colorRampPalette(c("dodgerblue4", "white", "#c0392b"))(paletteLength)
myBreaks <- c(seq(min(celltype220_data_H3K9ac), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(celltype220_data_H3K9ac)/paletteLength, max(celltype220_data_H3K9ac), length.out=floor(paletteLength/2)))
png("H3K9ac_z.png", height = 32, width  = 18.5, units="cm",res = 300)
pheatmap(
  mat               = celltype220_data_H3K9ac,
  color             = my.color,
  breaks           = myBreaks,
  cluster_rows=T,
  cluster_cols=F,
  annotation_row = group_H3K9ac,
  annotation_colors = annotation_colors,
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







##################test the difference between the 220 celltypes
baseline.diff <- TestDiff(baseline1.lua$Enrichment,baseline1.lua$Enrichment_std_error,baseline2.TN$Enrichment,baseline2.TN$Enrichment_std_error)
idx <- which(celltypes.diff<=0.05/53)
print(idx)
