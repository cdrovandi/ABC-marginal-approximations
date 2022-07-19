
# this code examines how the mean and variance of the median summary statistic is affected by the g-and-k parameters

library(gk)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(grid)

source("auxiliary_functions.R")

extra_args = list(nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10))
m = 50;

theta = matrix(0,nrow=1000,ncol=4)
summary = matrix(0,nrow=1000,ncol=m)
for (i in 1:1000){
  theta[i,] = prior_sim(extra_args)
  for (k in 1:m){
    x = rgk(extra_args$nobs,theta[i,1],theta[i,2],theta[i,3],theta[i,4])
    summary[i,k] = median(x)
  }
}

mean_summary = rowMeans(summary)

plot(theta[,1], mean_summary)
plot(theta[,2], mean_summary)
plot(theta[,3], mean_summary)
plot(theta[,4], mean_summary)

std_summary = rowSds(summary)

plot(theta[,1], std_summary)
plot(theta[,2], std_summary)
plot(theta[,3], std_summary)
plot(theta[,4], std_summary)


sims = data.frame(cbind(theta, mean_summary, std_summary))
colnames(sims) = c("a","b","g","k","mean_S1","sd_S1")

pa1 = ggplot(data=sims, aes(x=a, y=mean_S1)) + 
  geom_point() + theme(axis.title=element_text(size=16))

pb1 = ggplot(data=sims, aes(x=b, y=mean_S1)) + 
  geom_point() + theme(axis.title=element_text(size=16))

pg1 = ggplot(data=sims, aes(x=g, y=mean_S1)) + 
  geom_point() + theme(axis.title=element_text(size=16))

pk1 = ggplot(data=sims, aes(x=k, y=mean_S1)) + 
  geom_point() + theme(axis.title=element_text(size=16))


pa2 = ggplot(data=sims, aes(x=a, y=sd_S1)) + 
  geom_point() + theme(axis.title=element_text(size=16))

pb2 = ggplot(data=sims, aes(x=b, y=sd_S1)) + 
  geom_point() + theme(axis.title=element_text(size=16))

pg2 = ggplot(data=sims, aes(x=g, y=sd_S1)) + 
  geom_point() + theme(axis.title=element_text(size=16))

pk2 = ggplot(data=sims, aes(x=k, y=sd_S1)) + 
  geom_point() + theme(axis.title=element_text(size=16))


x11()
grid.arrange(pa1, pb1, pg1, pk1, pa2, pb2, pg2, pk2, nrow = 2)



