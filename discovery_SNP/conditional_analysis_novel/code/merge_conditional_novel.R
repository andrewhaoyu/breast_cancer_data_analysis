setwd('/data/zhangh24/breast_cancer_data_analysis/')
p.value.twostage <- rep(0,11)
p.value.standard <- rep(0,11)
for(i1 in 1:11){
  load(paste0("./discovery_SNP/conditional_analysis_novel/novel_conditional_reuslt",i1,".Rdata"))
  p.value.twostage[i1] <- result[[1]]
  p.value.standard[i1] <- result[[2]]
}
names(p.value.twostage) <- c("")




###########second time check
setwd('/data/zhangh24/breast_cancer_data_analysis/')
p.value.twostage <- rep(0,11)
p.value.standard <- rep(0,11)
for(i1 in 4:11){
  load(paste0("./discovery_SNP/conditional_analysis_novel/result/novel_conditional_result",i1,".Rdata"))
  p.value.twostage[i1] <- result[[1]]
  #p.value.standard[i1] <- result[[2]]
}
names(p.value.twostage) <- c("rs7760611",
                             "rs1027113",
                             "rs1061657",
                             "rs6697258:120485335:C:A",
                             "rs56826596:45374890:G:A",
                             "rs139331653:45939294:G:A",
                             "rs16988381:30592808:C:A",
                             "rs16988381:30592808:C:A",
                             "rs34044188:65257363:C:T",
                             "1:145126177:G:A",
                             "rs17743054:42900892:T:C")
