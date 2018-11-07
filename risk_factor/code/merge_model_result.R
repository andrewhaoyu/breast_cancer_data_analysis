setwd('/dcl01/chatterj/data/hzhang1/breast_cancer_data_analysis')
library(CompQuadForm)
GlobalTest <- function(theta,covar){
  GTA.stat <- t(theta)%*%solve(covar)%*%theta
  df <- length(theta)
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  places = 3
  power.number <- floor(-log10(p.value.GTA))+places
  ###format the output with three digits in total
  p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)
  
  # GTA.mixed <- theta%*%t(theta)
  # 
  # lamda <- eigen(covar)$values
  # 
  # result <- davies(GTA.mixed,lamda,lim = 2000000,acc=1e-9)
  # p.value.GTA.mixed <- result[[3]]
  
  
  return(c(p.value.GTA))
}
GlobalTestCate <- function(model,start,second.num,
                           covar.num){
  theta <- model[[1]][start+c(1:(second.num*covar.num))]
  covar <-   model[[2]][start+c(1:(second.num*covar.num)),
                        start+c(1:(second.num*covar.num))]
  result <- NULL
  for(i in 2:second.num){
    theta_temp<- theta[i+c(0:(covar.num-1))*second.num]
    covar_temp <- covar[i+c(0:(covar.num-1))*second.num,
                        i+c(0:(covar.num-1))*second.num]
    
    result <- rbind(result,GlobalTest(theta_temp,covar_temp))
  }
  if(second.num==5){
    row.names(result) <- c("ER","PR",
                           "HER2","grade")
    
  }else if(second.num==6){
    row.names(result) <- c("ER","PR",
                           "HER2","grade","TN interaction") 
  }
  return(result)  
} 
  i1= 1
  load(paste0("./risk_factor/result/model",i1,".Rdata"))
  start <- 23
  second.num <- 5
  covar.num <- 3
  ###########age first birth
  result1 <- GlobalTestCate(model,start,second.num,covar.num)
  #########breast_cat
  start <- start+second.num*(covar.num+1)
  second.num <- 5
  covar.num <- 4
  result2 <- GlobalTestCate(model,start,second.num,covar.num)
  
  #########parity
  start <- start+second.num*(covar.num+1)
  second.num <- 5
  covar.num <- 4
  result3 <- GlobalTestCate(model,start,second.num,covar.num)
  result <- t(cbind(result1,result2,result3))
  row.names(result) <- c("agefftp_cat",
                         "breastmos_cat",
                         "parity_cat")
library(xlsx)
  # write.xlsx(model[[4]]
  #            ,file = "risk_factor_result_110118.xlsx",
  #            sheetName = "pb_tr_study_continous",
  #            append = T)
  # write.xlsx(model[[5]]
  #            ,file = "risk_factor_result_110118.xlsx",
  #            sheetName = "pb_test_study_continous_tr",
  #            append=T)
  write.xlsx(model[[4]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "pb_additive",
             append=T) 
  write.xlsx(model[[5]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "pb_additive_test",
             append=T) 
  write.xlsx(result,
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "pb_additive_second_test",
             append=T) 
  
  
  i1= 2
  load(paste0("./risk_factor/result/model",i1,".Rdata"))
  start <- 23
  second.num <- 6
  covar.num <- 3
  ###########age first birth
  result1 <- GlobalTestCate(model,start,second.num,covar.num)
  #########breast_cat
  start <- start+second.num*(covar.num+1)
  second.num <- 6
  covar.num <- 4
  result2 <- GlobalTestCate(model,start,second.num,covar.num)
  
  #########parity
  start <- start+second.num*(covar.num+1)
  second.num <- 6
  covar.num <- 4
  result3 <- GlobalTestCate(model,start,second.num,covar.num)
  result <- t(cbind(result1,result2,result3))
  row.names(result) <- c("agefftp_cat",
                         "breastmos_cat",
                         "parity_cat")
  write.xlsx(model[[4]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "pb_TN",
             append=T) 
  write.xlsx(model[[5]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "pb_TN_test",
             append=T) 
  write.xlsx(result,
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "pb_TN_second_test",
             append=T) 
  i1= 3
  load(paste0("./risk_factor/result/model",i1,".Rdata"))
  write.xlsx(model[[4]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "pb_intrinsic",
             append=T) 
  write.xlsx(model[[5]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "pb_intrinsic_test",
             append=T) 
  
  
  
  i1= 4
  load(paste0("./risk_factor/result/model",i1,".Rdata"))
  start <- 23
  second.num <- 5
  covar.num <- 3
  ###########age first birth
  result1 <- GlobalTestCate(model,start,second.num,covar.num)
  #########breast_cat
  start <- start+second.num*(covar.num+1)
  second.num <- 5
  covar.num <- 4
  result2 <- GlobalTestCate(model,start,second.num,covar.num)
  
  #########parity
  start <- start+second.num*(covar.num+1)
  second.num <- 5
  covar.num <- 4
  result3 <- GlobalTestCate(model,start,second.num,covar.num)
  result <- t(cbind(result1,result2,result3))
  row.names(result) <- c("agefftp_cat",
                         "breastmos_cat",
                         "parity_cat")
 
  write.xlsx(model[[4]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "npb_additive",
             append=T) 
  write.xlsx(model[[5]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "npb_additive_test",
             append=T) 
  write.xlsx(result,
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "npb_additive_second_test",
             append=T) 
  
  
  i1= 5
  load(paste0("./risk_factor/result/model",i1,".Rdata"))
  start <- 23
  second.num <- 6
  covar.num <- 3
  ###########age first birth
  result1 <- GlobalTestCate(model,start,second.num,covar.num)
  #########breast_cat
  start <- start+second.num*(covar.num+1)
  second.num <- 6
  covar.num <- 4
  result2 <- GlobalTestCate(model,start,second.num,covar.num)
  
  #########parity
  start <- start+second.num*(covar.num+1)
  second.num <- 6
  covar.num <- 4
  result3 <- GlobalTestCate(model,start,second.num,covar.num)
  result <- t(cbind(result1,result2,result3))
  row.names(result) <- c("agefftp_cat",
                         "breastmos_cat",
                         "parity_cat")
  write.xlsx(model[[4]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "npb_TN",
             append=T) 
  write.xlsx(model[[5]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "npb_TN_test",
             append=T) 
  write.xlsx(result,
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "npb_TN_second_test",
             append=T) 
  i1= 6
  load(paste0("./risk_factor/result/model",i1,".Rdata"))
  write.xlsx(model[[4]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "npb_intrinsic",
             append=T) 
  write.xlsx(model[[5]],
             file = "./risk_factor/result/risk_factor_result_110518.xlsx",
             sheetName = "npb_intrinsic",
             append=T) 
  
  
  
  
  
  
  
  
  
  
  
  
  # rbind(model[[5]]
  #       )
  # 