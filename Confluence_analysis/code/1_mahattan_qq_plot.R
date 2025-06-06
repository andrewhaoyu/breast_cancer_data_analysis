source("/data/zhangh24//breast_cancer_data_analysis/Theme_publication_final/theme_publication_final.R")
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
# Load Regenie output
data <- fread("/data/DCEG_Confluence/JW/nextflow_GSA/results/Outcome.regenie.gz")

# Clean and format columns for plotting
data <- data %>%
  rename(
    CHR = CHROM,
    BP = GENPOS,
    rsid = ID,
    FREQ_A1 = A1FREQ
  ) %>%
  mutate(
    CHR = as.integer(CHR),
    BP = as.integer(BP),
    FREQ_A1 = as.numeric(FREQ_A1),
    INFO = as.numeric(INFO),
    P = as.numeric(10^(-LOG10P)),
    P = ifelse(P == 0, 1E-300, P),
    N = as.integer(N),
    MAF = ifelse(FREQ_A1 <= 0.5, FREQ_A1, 1 - FREQ_A1),
    MAC = MAF*N
  ) %>%
  filter(
    !is.na(CHR) & !is.na(BP) & !is.na(FREQ_A1) & !is.na(P) & !is.na(INFO),
    INFO > 0.2,
    MAC >= 30
  ) %>%
  select(rsid, CHR, BP, FREQ_A1, MAF, INFO, P, N, GENE_NAME)


dat = data %>%
  mutate(MAF = ifelse(FREQ_A1 <= 0.5, FREQ_A1, 1 - FREQ_A1)) %>%
  select(rsid, CHR, BP, P, MAF) %>%
  rename(SNP = rsid)

library(readr)


x = dat$P
z = qnorm(x / 2)
lambda = round(median(z^2) / qchisq(0.5,1), 3)
N.effect = median(data$N)
lambda_1000 = round(1+1000*(lambda-1)/N.effect  ,3)



convert.qval.pval = function(qvalues) {
  # you need to know the estimate of pi0 used to create the q-value
  # that's the maximum q-value (or very, very close to it)
  pi0 = max(qvalues)
  # compute m0, the estimated number of true nulls
  m0 = length(qvalues) * pi0
  # then you multiply each q-value by the proportion of true nulls
  # expected to be under it (the inverse of how you get there from
  # the p-value):
  return(qvalues * rank(qvalues) / m0)
}
p.pwas <- 5E-08
#q.pwas <- convert.qval.pval(c(tmp, 0.05))[(nrow(dat)+1)]

# Compute cumulative base pair position
nCHR <- length(unique(data$CHR))
offset <- 0  # space between chromosomes
s <- 0
nbp <- c()

for (i in sort(unique(data$CHR))) {
  max_bp <- max(data$BP[data$CHR == i])
  nbp[i] <- max_bp
  data$BPcum[data$CHR == i] <- data$BP[data$CHR == i] + s
  s <- s + max_bp + offset  # <-- add spacing here
}


# Prepare axis labels
axis.set <- data %>%
  group_by(CHR) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

# Define plotting limits
ylim <- abs(floor(log10(min(data$P)))) + 2

# Genome-wide significance line
sigline <- data.frame(sig = -log10(p.pwas), val = paste0("P=", signif(p.pwas, 2)))


# Identify lead SNPs for annotation (1 per ~500kb per chromosome)
top_hits <- data %>%
  filter(P < 1e-10) %>%
  arrange(CHR, BP) %>%
  group_by(CHR) %>%
  mutate(index = floor(BP / 3e6)) %>%  # Bin by 3 Mb
  group_by(CHR, index) %>%
  slice_min(order_by = P, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(logP = -log10(P)) %>%
  filter(!is.na(GENE_NAME)) %>%
  filter(!grepl("^AL|^RP|^LOC|^AC", GENE_NAME)) %>%  # remove generic names
  select(CHR, BPcum, logP, GENE_NAME)

top_hits <- top_hits %>%
  filter(!grepl("^AL|^RP|^LOC", GENE_NAME))  # optional cleanup


# Manhattan plot
manhplot <- ggplot(data, aes(x = BPcum, y = -log10(P), color = as.factor(CHR), size = -log10(P))) +
  geom_point(alpha = 0.8, size = 0.8) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#08306b", "#4292c6"), nCHR)) +
  scale_size_continuous(range = c(0.5, 3)) +
  geom_hline(data = sigline, aes(yintercept = sig), color = "red", linetype = "dashed") +
  geom_text_repel(
    data = top_hits,
    aes(x = BPcum, y = logP, label = GENE_NAME),
    size = 2.3,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.2,
    max.overlaps = Inf,
    nudge_y = 2,
    force = 2,
    min.segment.length = 0,
    seed = 42,
    show.legend = FALSE
  )+
  # geom_text(data = top_hits, aes(x = BPcum, y = logP + 1.5, label = GENE_NAME),
  #           size = 2.5, vjust = 0, check_overlap = TRUE, fontface = "bold") +
  guides(color = "none") +
  labs(x = NULL, y = "-log10(p)") +
  theme_Publication() +
  theme(
    legend.position = "top",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5),
    plot.title = element_text(size = 12, face = "bold")
  )


outpath <-"/data/zhangh24/breast_cancer_data_analysis/Confluence_analysis/result/"  


ggsave(filename=paste0("man_gsa.png"),
       plot=manhplot, device="png",
       path=outpath,
       width=9, height=4, units="in", dpi=300)


# Manhattan plot
manhplot_no_anno <- ggplot(data, aes(x = BPcum, y = -log10(P), color = as.factor(CHR), size = -log10(P))) +
  geom_point(alpha = 0.8, size = 0.8) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#08306b", "#4292c6"), nCHR)) +
  scale_size_continuous(range = c(0.5, 3)) +
  geom_hline(data = sigline, aes(yintercept = sig), color = "red", linetype = "dashed") +
  # geom_text(data = top_hits, aes(x = BPcum, y = logP + 1.5, label = GENE_NAME),
  #           size = 2.5, vjust = 0, check_overlap = TRUE, fontface = "bold") +
  guides(color = "none") +
  labs(x = NULL, y = "-log10(p)") +
  theme_Publication() +
  theme(
    legend.position = "top",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5),
    plot.title = element_text(size = 12, face = "bold")
  )


outpath <-"/data/zhangh24/breast_cancer_data_analysis/Confluence_analysis/result/"  


ggsave(filename=paste0("man_gsa_no_anno.png"),
       plot=manhplot_no_anno, device="png",
       path=outpath,
       width=9, height=4, units="in", dpi=300)

#qq(dat$P)










library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")



# option_list <- list(
#   make_option("--input", type="character", default="",
#               help="Input file, tab delimited; required columns: 'MAF' and 'PVALUE'"),
#   make_option("--prefix", type="character", default="",
#               help="Prefix of output files"),
#   make_option("--top.size", type="numeric", default=0.125,
#               help="top size = proportion of total length y axis [default=0.125]"),
#   make_option("--break.top", type="numeric", default=15,
#               help="set axis break at -log10(P) [default=15]"),
#   make_option("--width", type="numeric", default=900,
#               help="Width QQ plot in pixel [default=900]"),
#   make_option("--height", type="numeric", default=900,
#               help="Height QQ plot in pixel [default=900]"),
#   make_option("--pointsize", type="numeric", default=16,
#               help="Point size of plots [default=16]"),
#   make_option("--maf", type="character", default="",
#               help="name of column with MAF [default='']"),
#   make_option("--af", type="character", default="",
#               help="name of column with AF [default='']"),
#   make_option("--pvalue", type="character", default="PVALUE",
#               help="name of column with p.value [default='PVALUE']"),
#   make_option("--log10p", type="logical", default=F,
#               help="Input p.value column with -log10(p.value) [default=F]"),
#   make_option("--maintitle", type="character", default="",
#               help="Plot title")
# )
# 
# parser <- OptionParser(usage="%prog [options]", option_list=option_list)
# # 
# args <- parse_args(parser, positional_arguments = 0)
# opt <- args$options


qqplotdata <- function(logpvector){
  o = sort(logpvector,decreasing=T)
  e = -log10(ppoints(length(o)))       
  qqdata <- data.frame(o,e)
  qqdata$o <- round(qqdata$o,3)
  qqdata$e <- round(qqdata$e,3)
  keepU <- which(!duplicated(qqdata))
  qqdata <- qqdata[keepU,]
  
  N <- length(logpvector) ## number of p-values
  ## create the confidence intervals
  qqdata$c975 <- NA
  qqdata$c025 <- NA
  
  ## the jth order statistic from a
  ## uniform(0,1) sample
  ## has a beta(j,n-j+1) distribution
  ## (Casella & Berger, 2002,
  ## 2nd edition, pg 230, Duxbury)
  
  for(i in 1:length(keepU)){
    j <- keepU[i]
    qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
    qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
  }
  return(qqdata)
}
yLine <- c(-log10(5E-8))
colLine <- c("red")
dat$log10P = -log10(dat$P)
gwas = as.data.frame(dat)
# Determine frequency bins and create variable for binned QQ plot

minMAF <- min(gwas$MAF)

freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
gwas$freqbin <- cut(gwas$MAF, freqbins,include.lowest=T)
freqtable <- table(gwas$freqbin)
freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
freqtable <- freqtable[freqtable > 0]

## Generate QQ plot data by frequency bin
fbin <- character(0)
fN <- integer(0)
fx <- numeric(0)
fy <- numeric(0)
fcol <- character(0)
legendcol <- character(0)
conf <- list()
allcols <- brewer.pal(4,"Set1")
ycol <- "log10P"
for(f in 1:length(freqtable)){
  fbin <- c(fbin,names(freqtable)[f])
  fsnps <- which(gwas$freqbin ==names(freqtable)[f])
  plotdata <- qqplotdata(gwas[[ycol]][fsnps])
  fN <- c(fN,freqtable[f])
  fx <- c(fx,plotdata$e)
  fy <- c(fy,plotdata$o)
  fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
  conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                          'y'=c(plotdata$c975,rev(plotdata$c025)))
  legendcol <- c(legendcol,allcols[f])
}
legendtext <- paste0("MAF=",fbin,"; N SNPs=",format(fN,big.mark=",",scientific=FALSE))
opt =  list(break.top = 15,
            top.size = 0.125)


png(filename = paste0(outpath,"/QQ_gsa.png"), width = 8, height = 8, units = "in",res=300)
xlim <- c(0,max(fx,na.rm=T))
ylim <- c(0,max(fy,na.rm=T))
maxY <- max(fy,na.rm=T)
print("okkkk2")
par(mar=c(5.1,5.1,4.1,1.1))
print("okkkk3")

lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
lab2 <- lab2[lab2 > max(lab1)]

# resulting range of top scale in bottom scale units
top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
top.data = max(lab2)-opt$break.top

# function to rescale the top part
rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
rescaled.y = rescale(fy[fy>opt$break.top])
plot(0,0,
     ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
     xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
     ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
     cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
     main=opt$maintitle,pch=19)

# Plot confidence intervals	
for(p in 1:length(conf)){
  polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
          col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
          border = NA)
}

# add points
points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)

# identify line & add axis break
lines(xlim,xlim,col="black",lty = 2)
axis(1,cex.axis=1.5,cex.lab=1.5)
par(las=1)
axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
par(las=0)
box()
par(las=0)
points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
par(las=1)
axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
abline(h=ifelse(yLine<opt$break.top,
                yLine,
                rescale(yLine)),
       col=colLine,lwd=1.5,lty=2)
legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
text(5,1,expression(paste(lambda[1000]," = ")),cex = 1.5)
text(5.7,1,paste(lambda_1000),cex = 1.5)

title(paste0("GSA QQ plot"))
dev.off()

# plot(1,1)
# text(1,0.8,expression(paste(lambda[1000]," = ",buquote(.(lambda_1000)))),cex = 1.5)

# lambda_vec = c(lambda,lambda_1000)
# save(lambda_vec,file = paste0("/data/zhangh24/multi_ethnic/result/GLGC/lambda_value/lambda_vec_",i1,"_",i2,".rdata"))

