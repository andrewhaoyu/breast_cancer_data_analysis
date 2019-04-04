#goal: demonstrate box r in BCAC
#https://cran.r-project.org/web/packages/boxr/vignettes/boxr.html
#627lww8un9twnoa8f9rjvldf7kb56q1m
#gSKdYKLd65aQpZGrq9x4QVUNnn5C8qqm
library(boxr)
box_auth()
box_dir()

data <- box_read(419308975889)
library(ggplot2)
ggplot(data)+geom_bar(aes(AgeDiag1))
