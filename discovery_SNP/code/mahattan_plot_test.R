#test mahattan plot package to see which one looks best
library(qqman)
manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P" )
library(data.table)
library(dplyr)
data <- as.data.frame(fread("./discovery_SNP/result/ResultsMeta_GWAS_iCOGs_Onco_filter_R2_MAF.txt"))

gwas_result <- data %>%
  filter(p.meta<=10^-4) %>% 
  select(c(SNP.Onco,chr.Onco,Position.Onco,p.meta))
dim(gwas_result)
colnames(gwas_result) <- c("SNP","CHR","BP","P")
idx <- which(gwas_result$CHR==1&gwas_result$BP==110222901)

min(gwas_result$P)
fine_mapping <- read.csv("./data/filter_regions_standard.csv",header= T)
idx_cut <- NULL
start <- fine_mapping$start
end <- fine_mapping$end
CHR <- fine_mapping$CHR

#fine_mapping
for(i in 1:nrow(fine_mapping)){
  print(i)
  chr_temp <- CHR[i]
  start_temp <- start[i]
  end_temp <- end[i]
  idx <- which(gwas_result$CHR==chr_temp&gwas_result$BP>=start_temp&gwas_result$BP<=end_temp)
  idx_cut <- c(idx_cut,idx)
}
############duplicate variables won't mater
idx_cut <- unique(idx_cut)
gwas_result_filter <- gwas_result[-idx_cut,]
dim(gwas_result_filter)
min(gwas_result_filter$P)

gwas_result_filter[idx,]
idx.sig <- which(gwas_result_filter$P<=5*10^-8)
length(idx.sig)
head(gwas_result_filter[idx.sig,])
gwas_sig <- gwas_result_filter[idx.sig,]
head(gwas_result_filter[idx.sig,])
gwas_sig[150:550,]
manhattan(gwas_result_filter, col= c("blue4", "orange3"),
          cex = 0.6,suggestiveline = F)

library("CMplot")
CMplot(gwas_result_filter,plot.type="c")






library(colortools)

#col <- complementary("skyblue")
mahattanplot <- function(gwasResults){
  library(ggplot2)
  col <- c("navy","skyblue3")
  don <- gwasResults %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
    
    # Filter SNP to make the plot lighter
    filter(-log10(P)>0.5)
  
  # Prepare X axis
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  # Prepare text description for each SNP:
 # don$text <- paste("SNP: ", don$SNP, "\nPosition: ", don$BP, "\nChromosome: ", don$CHR, "\nLOD score:", -log10(don$P) %>% round(2), "\nWhat else do you wanna know", sep="")
  
  # Make the plot
  p <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=1, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    ylim(0,9) +
    
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
return(p)
  
}


