#goal: demonstrate box r in BCAC
#https://cran.r-project.org/web/packages/boxr/vignettes/boxr.htm?l
#627lww8un9twnoa8f9rjvldf7kb56q1m？
#gSKdYKLd？？ry(boxr)
box_auth()
#box_dir()

data <- box_read('419308975889')
data[1:10,1:10]
library(ggplot2)
ggplot(data)+geom_bar(aes(AgeDiag1))
