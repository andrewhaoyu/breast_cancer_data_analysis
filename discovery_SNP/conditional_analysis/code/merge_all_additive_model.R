result.com <- NULL








region <- 207
con.show <- rep(0,region)
same.all <- rep(0,region)
one.same.all <- rep(0,region)
each.other.same.all <- rep(0,region)



for(i1 in 1:207){
  print(i1)
  load(paste0("/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/result.all",i1,".Rdata"))
  if(is.null(result.all)==1){
    con.show[i1] <- 0
  }else{
    con.show[i1] <- 1
    n <- nrow(result.all)
    top.sig <- ifelse(result.all[1,c(10,12,14,16)]<=0.05,1,0)
    con.sig <- ifelse(result.all[2:n,c(10,12,14,16)]<=0.05,1,0)
    n.con <- nrow(con.sig)
    same <- 0
    one.same <- 0
    each.other.same <- 0
    if(n.con==1){
      each.other.same <- NA 
    }else{
      if(any(colSums(con.sig)>2)){
        each.other.same = 1
    }
    for(i in 1:n.con){
      if(all.equal(as.vector(con.sig[i,]),as.vector(top.sig))==1){
        same = same + 1
      }
      idx <- which(top.sig==1)
      if(any(con.sig[i,idx]==1)){
        one.same=one.same+1
      }
    }
    
    
    }
    same.all[i1] <- same
    one.same.all[i1] <- one.same
    each.other.same.all[i1] <- each.other.same
    
  }
}


write.csv(result.com,file="/spin1/users/zhangh24/breast_cancer_data_analysis/discovery_SNP/conditional_analysis/result/condition.all.csv")
colnames(result.com)[19]="known_flag"
table(result.com$mark)


