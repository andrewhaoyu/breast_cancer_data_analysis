#goal: merge dynamic prediction
#original coeficients
i3 = 1
pop.ind = 2
i1 = 1
n.s <- 100
r2.vad.mat <- matrix(0,n.s,4)
n.snp.mat <- matrix(0,n.s,4)
prop.mat <- matrix(0,n.s,4)
alpha.mat <- matrix(0,n.s,4)
p.thr.mat <- matrix(0,n.s,4)

#i3 is method i3=1 represent orignal coef
#i3 =2 represent normal prior coef
#i1 index of simulation
#pop.ind 2 is African
#pop.ind 3 is LAT
temp <- 1
for(i3 in 1:2){
  for(pop.ind in 2:3){
    
    for(i1 in 1:100){
      load(paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3))
      idx <- which.max(result[,3])
      p.thr.mat[i1,temp] <- result[idx,1]
      alpha.mat[i1,temp] <- result[idx,2]
      r2.vad.mat[i1,temp] <- result[idx,4]
      n.snp.mat[i1,temp] <- result[idx,5]
      prop.mat[i1,temp] <- result[idx,6]
      
    } 
    temp <- temp+1 
     }
 
}
r2.vad.mat
colMeans(r2.vad.mat)
colMeans(n.snp.mat)
colMeans(prop.mat)


#calibrated coefficents
i3 = 2



if(i3==1){
  result <- LDPDy(y_all,
                  beta.train,
                  p.train,
                  p.thr,
                  pop.ind,
                  beta.ref,
                  p.ref,
                  alpha)
  save(result,file=paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3))  
  
}else{
  result <- LDPDyW(y_all,
                   beta.train,
                   sd.train,
                   p.train,
                   p.thr,
                   pop.ind,
                   beta.ref,
                   sd.ref,
                   p.ref,
                   alpha)
  save(result,file=paste0("./multi_ethnic/result/Dy_result_",i1,"_",pop.ind,"_",i3))  
}            

