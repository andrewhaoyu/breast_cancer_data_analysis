#generate matrix illustration for Figure 1
library(reshape2)
library(ggplot2)
library(RColorBrewer)
source("./generate_plot/code/theme_publication.R")
nrow = 8
ncol = 10
set.seed(6)
A <- matrix(rbinom(nrow*ncol,2,0.33), 
            byrow = TRUE, nrow = nrow, ncol = ncol)
colnames(A) <- c(1:ncol)
rownames(A) <- c(1:nrow)
longData <- melt(A)
fill_color = brewer.pal(9, "Greens")[c(2,5,8)]


p1 = ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = as.character(value)), color = "black",size = 0.7) +  # Add black borders
  scale_fill_manual(values = fill_color)+
  labs(x = NULL, y = NULL, title = "Matrix") +  # Remove axis labels
  theme_Publication() +
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",  # Remove legend
    plot.title = element_blank(),
    axis.line = element_blank()  # Remove axis lines
  )
png(file = "./generate_plot/result/genotype_matrix.png",
    width = 3, height = 2, units = "in", res = 300)
print(p1)
dev.off()




nrow = 8
ncol = 1
A <- matrix(c(c(1,1,1,1),c(0,0,0,0)), 
            byrow = TRUE, nrow = nrow, ncol = ncol)
colnames(A) <- c(1:ncol)
rownames(A) <- c(1:nrow)
longData <- melt(A)
fill_color = brewer.pal(9, "Blues")[c(8,1)]


p2 = ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = as.character(value)), color = "black",size = 0.7) +  # Add black borders
  scale_fill_manual(values = fill_color)+
  labs(x = NULL, y = NULL, title = "Matrix") +  # Remove axis labels
  theme_Publication() +
  theme(
    plot.margin = margin(0, 0, 0, 0, "in"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",  # Remove legend
    plot.title = element_blank(),
    axis.line = element_blank()  # Remove axis lines
  )
png(file = "./generate_plot/result/phenotype_matrix.png",
    width = 0.375, height = 2, units = "in", res = 300)
print(p2)
dev.off()






nrow = 8
ncol = 4
A <- matrix(c(0,1,-9,1,888,888,888,888), 
            byrow = FALSE, nrow = nrow, ncol = ncol)
colnames(A) <- c(1:ncol)
rownames(A) <- c(1:nrow)
longData <- melt(A)
fill_color = brewer.pal(9, "Blues")[c(8,1)]


p2 = ggplot(longData, aes(x = Var2, y = Var1)) + 
  geom_tile(aes(fill = as.character(value)), color = "black",size = 0.7) +  # Add black borders
  scale_fill_manual(values = fill_color)+
  labs(x = NULL, y = NULL, title = "Matrix") +  # Remove axis labels
  theme_Publication() +
  theme(
    plot.margin = margin(0, 0, 0, 0, "in"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",  # Remove legend
    plot.title = element_blank(),
    axis.line = element_blank()  # Remove axis lines
  )
png(file = "./generate_plot/result/phenotype_matrix.png",
    width = 0.375, height = 2, units = "in", res = 300)
print(p2)
dev.off()
