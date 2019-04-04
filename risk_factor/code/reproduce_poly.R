#goal:reproduce jenny's group polytomous model reuslts, missing indicator methods

library(data.table)
library(bc2)
#data <- fread("./data/dataset_montse_20180522.txt")
setwd('/Users/zhangh24/GoogleDrive/breast_cancer_data_analysis')
data <- fread("./data/Dataset_Montse_20190322.txt")
#data <- fread("./data/Dataset_Montse_2018-10-10.txt")




##reproduce Jenny's group polytomous model
############collapes breast mos cat 0,1
############collapes age fftp cat 0,1
##############we only focus on the invasive nbreast cancer cases
data$status[data$status==2|data$status==3] <- 1
data$breastmos_cat[data$breastmos_cat==0] <- 1
data$agefftp_cat[data$agefftp_cat==0] <- 1
data$lastchildage_cat[data$lastchildage_cat==0] <- 1


##############only focus on population based study


library(dplyr)
data1 = data %>% filter(design_cat==0)
###############polytomous model 
library(nnet)
model1 <- multinom(molgroup~as.factor(agemenarche_cat)+
                     as.factor(parity_cat)+
                     as.factor(mensagelast_cat)+
                     as.factor(agefftp_cat)+
                     as.factor(breastmos_cat)+
                    # as.factor(lastchildage_cat)+
                     +study+refage,data=data1, maxit= 500)
coef.1 <- coef(model1)
covar.1 <- vcov(model1)
result <- list(coef.1,covar.1)
save(result,file = paste0("./risk_factor/result/poly/poly_result.Rdata"))
#interested in first 22 variable
n.in <- 22
colnames(coef.1)[1:n.in] <- c("Intercept",
                            paste0("agemenarche_cat",c(1,2,3,9)),
                            paste0("parity_cat",c(1,2,3,4,9)),
                            paste0("mensagelast_cat",c(1,2,9)),
                            paste0("agefftp_cat",c(2,3,4,9)),
                            paste0("breastmos_cat",c(2,3,4,5,9)))
                           # paste0("lastchildage_cat",c(2,3,4,9)))



#reorganize the covariance matrix to match coef data.frame
#columns are the variables
n.var <- ncol(coef.1)
#rows are the 5 intrinsic subtypes plus 9 
n.sub <- 6
var.1 <- coef.1
for(i in 1:n.var){
  var.1[,i] <- diag(covar.1[c(i+n.var*(0:(n.sub-1))),c(i+n.var*(0:(n.sub-1)))])
}

result <- GenerateP(coef.1,var.1)
p <- result[[1]]
low.95 <- result[[2]]
high.95 <- result[[3]]

coef.1.n <- coef.1[,2:n.in]
var.1.n <- var.1[,2:n.in]
p.n <- p[,2:n.in]
low.95.n <- low.95[,2:n.in]
high.95.n <- high.95[,2:n.in]

coef.l <- as.vector(t(coef.1.n))
var.l <- as.vector(t(var.1.n))
p.l <- as.vector(t(p.n))
low.95.l <- as.vector(t(low.95.n))
high.95.l <- as.vector(t(high.95.n))
variable.l <- rep(colnames(coef.1.n),n.sub)
subtypes.l <- c(rep("Luminal A-like",ncol(coef.1.n)),
                  rep("Luminal B,HER2-negative-like",ncol(coef.1.n)),
                  rep("Luminal B-like",ncol(coef.1.n)),
                  rep("HER2 enriched-like",ncol(coef.1.n)),
                  rep("TN",ncol(coef.1.n)),
                  rep("missing",ncol(coef.1.n)))
#create the variable category for colour
variable.cat <- substr(variable.l,1,nchar(variable.l)-1)

new.data <- data.frame(variable.l,subtypes.l,
                       coef.l,
                       var.l,
                       p.l,
                       low.95.l,
                       high.95.l,
                       exp(coef.l),
                       exp(low.95.l),
                       exp(high.95.l),
                       variable.cat,
                       stringsAsFactors = F)
colnames(new.data) <- c("variable",
                        "subtypes",
                        "coef",
                        "variance",
                        "p-value",
                        "low95.coef",
                        "high95.coef",
                        "OR",
                        "ORlow95",
                        "ORhigh95",
                        "variable_cat")
###########take out all the missing data for comparasion
new.data.c <- new.data[-grep("9",variable.l),]
subtypes <- c("Luminal A-like","Luminal B,HER2-negative-like",
              "Luminal B-like",
              "HER2 enriched-like",
              "TN")
for(i in 1:length(subtypes)){
  subtype.temp = subtypes[i]
  png(paste0("./risk_factor/result/poly/poly_result_",subtype.temp,".png"),height=20,width = 15,res=300,units="cm")
  new.data.plot <- new.data.c %>% filter(subtypes==subtype.temp)
  print(
    ggplot(new.data.plot,aes(x=variable,y=OR,ymin=ORlow95,ymax=ORhigh95,colour=variable_cat))+
    geom_pointrange()+
    #scale_colour_manual(values=c("#386cb0","#fdb462"))+
    #geom_line(yintercept=1,lty=2)+
    coord_flip()+
    theme_bw()+
    geom_hline(yintercept=1,size=1,lty=2)+
    theme(legend.position="none")+
    ylab("Odds ratio")+
    #scale_y_continuous(breaks=c(0,1,2.5,5,8))+
    #xlab("SNP")+
    # facet_grid(.~method)+
    ggtitle(paste0("Forest plot for odds ratio of ",subtype.temp)
    )+
    theme(plot.title = element_text(hjust=0.5,face="bold"),
          axis.text=element_text(face="bold"),
          axis.title.x = element_text(face="bold"),
          axis.title.y = element_text(face="bold"))
  )
  dev.off()
  
}



GenerateP <- function(coef,sigma){
  z <- coef/sqrt(sigma)
  p <- 2*pnorm(-abs(z),lower.tail = T)
  low.95 <- coef-1.96*sqrt(sigma)
  high.95 <- coef+1.96*sqrt(sigma)
  return(list(p,low.95,high.95))
}













#delete <- na.action(model1)
#print(delete)


######puting all of the missing data as NA for mice package to run
agemenarche_cat <- data1$agemenarche_cat
idx <- which(agemenarche_cat==9)
agemenarche_cat[idx] <- NA
parity_cat <- data1$parity_cat
idx <- which(parity_cat==9)
parity_cat[idx] <- NA
mensagelast_cat <- data1$mensagelast_cat
idx <- which(mensagelast_cat==9)
mensagelast_cat[idx] <- NA
agefftp_cat <- data1$agefftp_cat
idx <- which(agefftp_cat==9)
agefftp_cat[idx] <- NA
breastmos_cat <- data1$breastmos_cat
idx <- which(breastmos_cat==9)
breastmos_cat[idx] <- NA
lastchildage_cat <- data1$lastchildage_cat
idx <- which(lastchildage_cat==9)
lastchildage_cat[idx] <- NA
study <- data1$study
refage <- data1$refage
all.covariates <- data.frame(as.factor(agemenarche_cat),
                               as.factor(parity_cat),
                               as.factor(mensagelast_cat),
                               as.factor(agefftp_cat),
                               as.factor(breastmos_cat),
                               as.factor(lastchildage_cat),
                               study,
                              refage)

# time1 = proc.time()
# imp <- mice(all.covariates,m=1,seed=1,print=FALSE)
# time = proc.time()-time1
# head(imp$imp$as.factor.agemenarche_cat.)
# all.covariates.c <- complete(imp,1)
# 
# x <- seq(0,1,0.01)
# y <- 1+sin(x)
# plot(x,y)
# data <- data.frame(x,y)
# library(ggplot2)
# ggplot(data)+geom_point(aes(x,y))
# M <- 10
# y_new <- 1+y/M
# ggplot(data)+geom_point(aes(x,y_new))
# 
# 













###########create the dummy variable for agemenarchecat
###########create the dummy variable for parity_cat
###########create the dummy variable for mensagelast_cat
###########create the dummy variable for breastmos_cat
###########create the dummy variable for agefftp_cat
###########create the dummy variable for lastchildage_mat
agemenarche_mat <- model.matrix(~as.factor(data1$agemenarche_cat)-1)[,-1]
colnames(agemenarche_mat) <- paste0("agemenarche_cat",c(1:3,9))
parity_mat <- model.matrix(~as.factor(data1$parity_cat)-1)[,-1]
colnames(parity_mat) <- paste0("parity_cat",c(1:4,9))
mensagelast_mat <- model.matrix(~as.factor(data1$mensagelast_cat)-1)[,-1]
colnames(mensagelast_mat) <- paste0("mensagelast_cat",c(1:2,9))
agefftp_mat <- model.matrix(~as.factor(data1$agefftp_cat)-1)[,-1]
colnames(agefftp_mat) <- paste0("agefftp",c(1:4,9))
breastmos_mat <- model.matrix(~as.factor(data1$breastmos_cat)-1)[,-1]
colnames(breastmos_mat) <- paste0("breastmos_cat",c(1:5,9))
lastchildage_mat <- model.matrix(~as.factor(data1$lastchildage_cat)-1)[,-1]
colnames(lastchildage_mat) <- paste0("lastchildage_cat",c(1:4,9))
refage <- data$refage





# idx.try <- which(data$design_cat==0&data$molgroup==2)
# data.new <- data[idx.try,]
# table(data.new$ER_status1,data.new$PR_status1,
#       data.new$HER2_status1,data.new$Grade1)

#############put the missing tumor characteristics as 888
idx.ER.mis <- which(data$status==1&is.na(data$ER_status1))
data$ER_status1[idx.ER.mis] <- 888
idx.PR.mis <- which(data$status==1&is.na(data$PR_status1))
data$PR_status1[idx.PR.mis] <- 888
idx.HER2.mis <- which(data$status==1&is.na(data$HER2_status1))
data$HER2_status1[idx.HER2.mis] <- 888
idx.grade.mis <- which(data$status==1&is.na(data$Grade1))
data$Grade1[idx.grade.mis] <- 888
#############put the subject BCAC-16687664 HER2 status as 1
idx <- which(data$HER2_status1==2)
data$HER2_status1[idx] <- 1
#############check the result
table(data$status,data$ER_status1)
table(data$status,data$PR_status1)
table(data$status,data$HER2_status1)
table(data$status,data$Grade1)
library(nnet)
############collapes breast mos cat 0,1
############collapes age fftp cat 0,1
data$breastmos_cat[data$breastmos_cat==0] <- 1
data$agefftp_cat[data$agefftp_cat==0] <- 1
###########create the phenotype file
y <- cbind(data$status,data$ER_status1,data$PR_status1,
           data$HER2_status1,data$Grade1)
colnames(y) <- c("casecontrol",
                 "ER",
                 "PR",
                 "HER2",
                 "Grade")

############don't adjust for study
############population based study
idx1 <- which(data$design_cat==0)
###########two-stage model based on population based study
model.1 <- TwoStageModel(y=y[idx1,],
                         additive = cbind(agefftp_mat,
                                         breast_mat,
                                         parity_mat,
                                         refage)[idx1,],
                         missingTumorIndicator = 888
                         )
write.xlsx(model.1[[4]]
  ,file = "risk_factor_result_110118.xlsx",
           sheetName = "population_based_second_stage")
write.xlsx(model.1[[5]]
           ,file = "risk_factor_result_110118.xlsx",
           sheetName = "population_based_second_stage")


############ population based study




# model.1 <- multinom(molgroup~as.factor(agefftp_cat)+as.factor(breastmos_cat)+as.factor(parity_cat)+study+refage,data=data1, maxit= 500)
# coef(model.1)




############non population based study
idx2 <- which(data$design_cat==1)
data2 <- data[idx2,]
model.2 <- multinom(molgroup~as.factor(breastmos_cat)+study+refage,data=data2, maxit= 500)
coef(model.2)




