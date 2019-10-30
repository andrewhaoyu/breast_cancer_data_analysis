#-------------------------------------------------------------------
# Update Date: 11/26/2018
# Create Date: 11/26/2018
# Goal: generate prognosis score for different subtypes
# Author: Haoyu Zhang
#-------------------------------------------------------------------
#---------------------------------------#---------------------------------------

library(tidyverse)

rm(list=ls())

# This requires a data frame with the following variables:

# Variables for which missing data not permitted

# age.start (age at diagnosis), 
# screen (mode of detection, 0=clinical, 1=screening),
# grade (pathology grade, 1-3)
# nodes (number of positive nodes, micrometastases=0.5)
# size (max size in mm), er (er status, 0=negative, 1=positive)
# generation (adjuvant chemotherapy, 0=none, 2=2nd gen, 3=3rd gen)
# horm (adjuvant hormone therapy, 0=none, 1=yes)
# traz (trastuzumab therapy, 0=none, 1=yes)
# bis (bisphosphoate therapy, 0=none, 1=yes)
# time (This is the follow-up time in years for patients who have survived and
#       the hypothetical follow-up time that would have been observed for patients
#       that died if they had survived to end of follow up.  It is important that 
#       this is NOT just time to death.)

# Variables for which missing values allowed

# her2 (her2 status, 0=negative, 1=positive, 9=missing)
# ki67 (ki67 status, 0=negative, 1=positive, 9=missing)

# Creating an example data set
age.start <- seq(25,75,10)
screen <- c(0, 1)
size <- c(5, 20, 30, 50)
grade <- c(1, 2, 3)
nodes <- c(0, 0.5, 3, 10)
er <- c(0, 1)
her2 <- c(0, 1)
ki67 <- c(0, 1)
generation <- c(0, 2, 3)
horm <- c(0, 1)
traz <- c(0, 1)
bis <- c(0, 1)

df <- as.tibble(expand.grid(age.start = age.start, 
                            grade = grade, 
                            nodes = nodes, 
                            size = size, 
                            er = er, 
                            her2 = her2, 
                            ki67 = ki67, 
                            screen = screen,
                            generation = generation,
                            horm = horm,
                            traz = traz,
                            bis  =bis)) 
time <- 10 
# Random number for selecting subset of dataset
inc <- runif(length(df$er))

df <- cbind(df, time, inc)
df <- filter(df, inc>0.9)
df$br <- NA
df$oth <- NA
df$all <- NA

# Count the number of rows in the dataset for looping through calculations
rows <- length(df$er)

# Carry out the PREDICT calculations for each row

for (n in 1:rows) {
  
  # Input case characteristics
  age.start  <- df$age.start[n]
  screen     <- df$screen[n]     # Clinically detected = 0, screen detected = 1
  size       <- df$size[n]   # Tumour size mm
  grade      <- df$grade[n]     # Tumour grade
  nodes      <- df$nodes[n]     # Number positive nodes. Nodal micrometastases = 0.5
  er         <- df$er[n]     # ER+ = 1, ER- = 0
  her2       <- df$her2[n]     # HER2+ = 1, HER2- = 0, missing = 9
  ki67       <- df$ki67[n]     # KI67+ = 1, KI67- = 0, missing = 9
  generation <- df$generation[n]     # Chemo generation 0, 2 or 3 only
  horm       <- df$horm[n]     # Hormone therapy Yes = 1, no = 0
  traz       <- df$traz[n]    # Trastuzumab therapy Yes = 1, no = 0
  bis        <- df$bis[n]     # Bisphosphonate therapy Yes = 1, no = 0
  
  # Grade variable for ER neg
  grade.val <- ifelse(er==1, grade, ifelse(grade==1, 0, 1)) 
  
  # Generate the coefficients
  
  age.mfp.1   <- ifelse(er==1, (age.start/10)^-2-.0287449295, age.start-56.3254902)
  age.beta.1  <- ifelse(er==1, 34.53642, 0.0089827)
  age.mfp.2   <- ifelse(er==1, (age.start/10)^-2*log(age.start/10)-.0510121013, 0)
  age.beta.2  <- ifelse(er==1, -34.20342, 0)
  size.mfp    <- ifelse(er==1, log(size/100)+1.545233938, (size/100)^.5-.5090456276)
  size.beta   <- ifelse(er==1, 0.7530729, 2.093446)
  nodes.mfp   <- ifelse(er==1,log((nodes+1)/10)+1.387566896,
                        log((nodes+1)/10)+1.086916249)
  nodes.beta  <- ifelse(er==1, 0.7060723, .6260541)
  grade.beta  <- ifelse(er==1, 0.746655, 1.129091)
  screen.beta <- ifelse(er==1, -0.22763366, 0)
  her2.beta   <- ifelse(her2==1, 0.2413,
                        ifelse(her2==0, -0.0762,0 ))
  ki67.beta   <- ifelse(ki67==1 & er==1, 0.14904,
                        ifelse(ki67==0 & er==1, -0.1133,0 ))
  
  # Calculate the other and breast mortality indicies
  
  # Other mortality prognostic index (mi)
  mi <- 0.0698252*((age.start/10)^2-34.23391957)
  
  # Breast cancer mortality prognostic index (pi)
  pi <- age.beta.1*age.mfp.1 + age.beta.2*age.mfp.2 + size.beta*size.mfp +
    nodes.beta*nodes.mfp + grade.beta*grade.val + screen.beta*screen + 
    her2.beta + ki67.beta
  
  # Treatment coefficients
  c     <- ifelse(generation == 0, 0, ifelse(generation == 2, -.248, -.446))
  h     <- ifelse(horm==1 & er==1, -0.3857, 0)
  t     <- ifelse(her2==1 & traz==1, -.3567, 0)
  b     <- ifelse(bis==1, -0.198, 0) # Only applicable to menopausal women.
  
  # Treatment combined
  rx <- h + c + t + b
  
  # Non breast cancer mortality
  # Generate cumulative baseline other mortality
  base.m.cum.oth <- exp(-6.052919 + (1.079863*log(df$time[n])) + (.3255321*df$time[n]^.5))
  
  # Generate cumulative survival non-breast mortality
  s.cum.oth <- exp(-exp(mi)*base.m.cum.oth)
  
  # Generate annual survival from cumulative survival
  m.cum.oth <- 1 - s.cum.oth 
  
  
  # Breast cancer specific mortality
  # Generate cumulative baseline breast mortality
  if (er==1) {
    base.m.cum.br <- exp(0.7424402 - 7.527762/df$time[n]^.5 - 1.812513*log(df$time[n])/df$time[n]^.5) 
  } else { base.m.cum.br <- exp(-1.156036 + 0.4707332/df$time[n]^2 - 3.51355/df$time[n])
  }
  
  # Calculate the cumulative breast cancer survival
  s.cum.br <- exp(-exp(pi+rx)*base.m.cum.br)
  m.cum.br <- 1 - s.cum.br
  
  # All cause mortality
  m.cum.all <- 1 - s.cum.oth*s.cum.br
  s.cum.all <- 100-100*m.cum.all
  
  
  # Proportion of all cause mortality that is breast cancer
  prop.br <- m.cum.br/(m.cum.br+m.cum.oth)
  prop.oth <- m.cum.oth/(m.cum.br+m.cum.oth)
  
  # Predicted cumulative breast specific mortality
  pred.m.br    <- prop.br*m.cum.all
  
  # Predicted cumulative non-breast cancer mortality
  pred.m.oth <- prop.oth*m.cum.all
  
  # Predicted cumulative all-cause mortality
  pred.all <- pred.m.br + pred.m.oth
  
  # Predicted breast cancer mortality
  df$br[n] <- pred.m.br
  
  # Predicted non-breast cancer mortality
  df$oth[n] <- pred.m.oth
  
  # Predicted all cause mortality
  df$all[n] <- pred.all
  
}



SubtypesTrans <- function(casecon,ER,PR,HER2,grade){
  n <- length(casecon)
  result <- rep("unknown",n)
  idx.con <- which(casecon==0)
  result[idx.con] <- "control"
  idx.LA <- which(HER2==0&(ER==1|PR==1)&(grade==1|grade==2))
  idx.LB <- which(HER2==1&(ER==1|PR==1))
  idx.LUBHER2 <- which(HER2==0&(ER==1|PR==1)&grade==3)
  idx.HER2 <- which(HER2==1&ER==0&PR==0)
  idx.TN <- which(HER2==0&ER==0&PR==0)
  #idx.mis <- which(HER2==888|ER==888|PR==888|Grade==888)
  result[idx.LA] <- "Luminal_A"
  result[idx.LB] <- "Luminal_B"
  result[idx.LUBHER2] <- "Luminal_B_HER2Neg"
  result[idx.HER2] <- "HER2Enriched"
  result[idx.TN] <- "TripleNeg"
  result <- factor(result,levels=c("control",
                                   "Luminal_A",
                                   "Luminal_B",
                                   "Luminal_B_HER2Neg",
                                   "HER2Enriched",
                                   "TripleNeg",
                                   "unknown"))
  
  return(result)
  
}
library(R.utils)
library(data.table)
setwd("/data/zhangh24/breast_cancer_data_analysis/")



data1 <- fread("./data/PRS_subtype_icgos_pheno_v10_euro.csv",header=T)
data1 <- as.data.frame(data1)
data2 <- fread("./data/PRS_subtype_Onco_euro_v10_08012018.csv",
               header=T)
data2 <- as.data.frame(data2)
library(xlsx)
################take out all the people with in-situ status
data1 = data1[data1$Behaviour1!=2,]
data2 = data2[data2$Behaviour1!=2,]
