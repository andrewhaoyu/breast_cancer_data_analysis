install_github("andrewhaoyu/TOP")
library(devtools)
library(TOP)
data(data, package="TOP")
model.fixed <- fixedEffect(data,"Behavior",
                     c("ER","PR"),"SNP",
                     covars.obj = ~PC1+PC2)
model[[1]]
modelInfo(data, "Behavior", c("ER", "PR"), covars.obj=~PC1 + PC2)
model.random <- randomEffect(data,"Behavior",
            c("ER","PR"),"SNP",
            covars.obj = ~PC1+PC2)


gfile <- system.file("sampleData", "chr1_1.impute.txt", package="TOP")
sfile <- system.file("sampleData", "subjects.impute.txt", package="TOP")
pfile <- system.file("sampleData", "pheno.txt", package="TOP")
# Define the list for the genotype data.
geno.list <- list(file=gfile, format="impute", subject.list=sfile)
# Only process the first 5 SNPs in the file
geno.list$start.vec <- 1
geno.list$stop.vec <- 6
# Define pheno.list
pheno.list <- list(file=pfile, delimiter="\t", id.var="ID")
# Define the variables in the model: Status ~ Age
pheno.list$response.var <- "Status"
pheno.list$main.vars <- "Age"
scan.score(geno.list, pheno.list)
 
