

# Package names
packages <- c("ggplot2", "gridExtra", "grid")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# load libraries
#library(ggplot2)
#library(gridExtra)
#library(grid)

source("auxiliary_functions.R")

# load the "observed" data
load(file = "data_gandk.RData")
obs_summ = summary_robust(y)

true_param = c(3,1,2,0.5)


##################### parameter 1 #####################

load("results_robust.RData")
theta_robust = 10/(1+exp(-results$theta))
theta_robust = data.frame(theta_robust)
theta_robust$method = "all summs"
colnames(theta_robust) = c("a", "b", "g", "k", "method")

load("results_robust_pilot.RData")
theta_robust_pilot = 10/(1+exp(-results$theta))
theta_robust_pilot = data.frame(theta_robust_pilot)
theta_robust_pilot$method = "pilot (all summs)"
colnames(theta_robust_pilot) = c("a", "b", "g", "k", "method")

load("results_robust1.RData")
theta_robust1 = 10/(1+exp(-results$theta))
theta_robust1 = data.frame(theta_robust1)
theta_robust1$method = "summ 1"
colnames(theta_robust1) = c("a", "b", "g", "k", "method")

load("results_robust1_continue.RData")
theta_robust1_continue = 10/(1+exp(-results$theta))
theta_robust1_continue = data.frame(theta_robust1_continue)
theta_robust1_continue$method = "pilot + summ 1"
colnames(theta_robust1_continue) = c("a", "b", "g", "k", "method")



theta = rbind(theta_robust,theta_robust1,theta_robust_pilot,theta_robust1_continue)

p1 = ggplot(data=theta) + 
  geom_density(aes(x=a, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position="none") +
  geom_point(aes(x=true_param[1], y=0), colour="black", size = 4) +
  xlim(c(2.75,3.25)) + theme(axis.title=element_text(size=16)) + labs(y= "density")


##################### parameter 2 #####################

load("results_robust.RData")
theta_robust = 10/(1+exp(-results$theta))
theta_robust = data.frame(theta_robust)
theta_robust$method = "all summs"
colnames(theta_robust) = c("a", "b", "g", "k", "method")

load("results_robust_pilot.RData")
theta_robust_pilot = 10/(1+exp(-results$theta))
theta_robust_pilot = data.frame(theta_robust_pilot)
theta_robust_pilot$method = "pilot (all summs)"
colnames(theta_robust_pilot) = c("a", "b", "g", "k", "method")

load("results_robust2.RData")
theta_robust2 = 10/(1+exp(-results$theta))
theta_robust2 = data.frame(theta_robust2)
theta_robust2$method = "summ 2"
colnames(theta_robust2) = c("a", "b", "g", "k", "method")

load("results_robust2_continue.RData")
theta_robust2_continue = 10/(1+exp(-results$theta))
theta_robust2_continue = data.frame(theta_robust2_continue)
theta_robust2_continue$method = "pilot + summ 2"
colnames(theta_robust2_continue) = c("a", "b", "g", "k", "method")



theta = rbind(theta_robust,theta_robust2,theta_robust_pilot,theta_robust2_continue)

p2 = ggplot(data=theta) + 
  geom_density(aes(x=b, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position="none") +
  geom_point(aes(x=true_param[2], y=0), colour="black", size = 4) +
  xlim(c(0.5,1.5)) + theme(axis.title=element_text(size=16)) + labs(y= "density")




##################### parameter 3 #####################



load("results_robust.RData")
theta_robust = 10/(1+exp(-results$theta))
theta_robust = data.frame(theta_robust)
theta_robust$method = "all summs"
colnames(theta_robust) = c("a", "b", "g", "k", "method")

load("results_robust_pilot.RData")
theta_robust_pilot = 10/(1+exp(-results$theta))
theta_robust_pilot = data.frame(theta_robust_pilot)
theta_robust_pilot$method = "pilot (all summs)"
colnames(theta_robust_pilot) = c("a", "b", "g", "k", "method")

load("results_robust3.RData")
theta_robust3 = 10/(1+exp(-results$theta))
theta_robust3 = data.frame(theta_robust3)
theta_robust3$method = "summ 3"
colnames(theta_robust3) = c("a", "b", "g", "k", "method")

load("results_robust3_continue.RData")
theta_robust3_continue = 10/(1+exp(-results$theta))
theta_robust3_continue = data.frame(theta_robust3_continue)
theta_robust3_continue$method = "pilot + summ 3"
colnames(theta_robust3_continue) = c("a", "b", "g", "k", "method")



theta = rbind(theta_robust,theta_robust3,theta_robust_pilot,theta_robust3_continue)

p3 = ggplot(data=theta) + 
  geom_density(aes(x=g, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position="none") +
  geom_point(aes(x=true_param[3], y=0), colour="black", size = 4) +
  xlim(c(1,3)) + theme(axis.title=element_text(size=16)) + labs(y= "density")


##################### parameter 4 #####################



load("results_robust.RData")
theta_robust = 10/(1+exp(-results$theta))
theta_robust = data.frame(theta_robust)
theta_robust$method = "all summs"
colnames(theta_robust) = c("a", "b", "g", "k", "method")

load("results_robust_pilot.RData")
theta_robust_pilot = 10/(1+exp(-results$theta))
theta_robust_pilot = data.frame(theta_robust_pilot)
theta_robust_pilot$method = "pilot (all summs)"
colnames(theta_robust_pilot) = c("a", "b", "g", "k", "method")

load("results_robust4.RData")
theta_robust4 = 10/(1+exp(-results$theta))
theta_robust4 = data.frame(theta_robust4)
theta_robust4$method = "summ"
colnames(theta_robust4) = c("a", "b", "g", "k", "method")

load("results_robust4_continue.RData")
theta_robust4_continue = 10/(1+exp(-results$theta))
theta_robust4_continue = data.frame(theta_robust4_continue)
theta_robust4_continue$method = "pilot + summ"
colnames(theta_robust4_continue) = c("a", "b", "g", "k", "method")



theta = rbind(theta_robust,theta_robust4,theta_robust_pilot,theta_robust4_continue)

p4 = ggplot(data=theta) + 
  geom_density(aes(x=k, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position = c(0.8, 0.8)) +
  theme(legend.title=element_blank()) +
  geom_point(aes(x=true_param[4], y=0), colour="black", size = 4) +
  xlim(c(0.25,0.75)) + theme(axis.title=element_text(size=16)) + labs(y= "density")



x11()
grid.arrange(p1, p2, p3, p4, nrow = 1)






################## s1 ##############


load("results_robust.RData")
ssx = results$ssx[,1]
ssx = data.frame(ssx)
ssx$method = "all summs"
colnames(ssx) = c("s1", "method")

load("results_robust_pilot.RData")
ssx_pilot = results$ssx[,1]
ssx_pilot = data.frame(ssx_pilot)
ssx_pilot$method = "pilot (all summs)"
colnames(ssx_pilot) = c("s1", "method")

load("results_robust1.RData")
ssx1 = results$ssx
ssx1 = data.frame(ssx1)
ssx1$method = "summ 1"
colnames(ssx1) = c("s1", "method")

load("results_robust1_continue.RData")
ssx1_continue = results$ssx
ssx1_continue = data.frame(ssx1_continue)
ssx1_continue$method = "pilot + summ 1"
colnames(ssx1_continue) = c("s1", "method")



ssx_results = rbind(ssx,ssx_pilot,ssx1,ssx1_continue)

p1 = ggplot(data=ssx_results) + 
  geom_density(aes(x=s1, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position="none") +
  geom_point(aes(x=obs_summ[1], y=0), colour="black", size = 4) +
  xlim(c(2.97,3.02)) + theme(axis.title=element_text(size=16)) + labs(y= "density")



################## s2 ##############


load("results_robust.RData")
ssx = results$ssx[,2]
ssx = data.frame(ssx)
ssx$method = "all summs"
colnames(ssx) = c("s2", "method")

load("results_robust_pilot.RData")
ssx_pilot = results$ssx[,2]
ssx_pilot = data.frame(ssx_pilot)
ssx_pilot$method = "pilot (all summs)"
colnames(ssx_pilot) = c("s2", "method")

load("results_robust2.RData")
ssx2 = results$ssx
ssx2 = data.frame(ssx2)
ssx2$method = "summ 2"
colnames(ssx2) = c("s2", "method")

load("results_robust2_continue.RData")
ssx2_continue = results$ssx
ssx2_continue = data.frame(ssx2_continue)
ssx2_continue$method = "pilot + summ 2"
colnames(ssx2_continue) = c("s2", "method")



ssx_results = rbind(ssx,ssx_pilot,ssx2,ssx2_continue)

p2 = ggplot(data=ssx_results) + 
  geom_density(aes(x=s2, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position="none") +
  geom_point(aes(x=obs_summ[2], y=0), colour="black", size = 4) +
  xlim(c(1.64,1.68)) + theme(axis.title=element_text(size=16)) + labs(y= "density")





################## s3 ##############


load("results_robust.RData")
ssx = results$ssx[,3]
ssx = data.frame(ssx)
ssx$method = "all summs"
colnames(ssx) = c("s3", "method")

load("results_robust_pilot.RData")
ssx_pilot = results$ssx[,3]
ssx_pilot = data.frame(ssx_pilot)
ssx_pilot$method = "pilot (all summs)"
colnames(ssx_pilot) = c("s3", "method")

load("results_robust3.RData")
ssx3 = results$ssx
ssx3 = data.frame(ssx3)
ssx3$method = "summ 3"
colnames(ssx3) = c("s3", "method")

load("results_robust3_continue.RData")
ssx3_continue = results$ssx
ssx3_continue = data.frame(ssx3_continue)
ssx3_continue$method = "pilot + summ 3"
colnames(ssx3_continue) = c("s3", "method")



ssx_results = rbind(ssx,ssx_pilot,ssx3,ssx3_continue)

p3 = ggplot(data=ssx_results) + 
  geom_density(aes(x=s3, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position="none") +
  geom_point(aes(x=obs_summ[3], y=0), colour="black", size = 4) +
  xlim(c(0.46,0.50)) + theme(axis.title=element_text(size=16)) + labs(y= "density")






################## s4 ##############


load("results_robust.RData")
ssx = results$ssx[,4]
ssx = data.frame(ssx)
ssx$method = "all summs"
colnames(ssx) = c("s4", "method")

load("results_robust_pilot.RData")
ssx_pilot = results$ssx[,4]
ssx_pilot = data.frame(ssx_pilot)
ssx_pilot$method = "pilot (all summs)"
colnames(ssx_pilot) = c("s4", "method")

load("results_robust4.RData")
ssx4 = results$ssx
ssx4 = data.frame(ssx4)
ssx4$method = "summ"
colnames(ssx4) = c("s4", "method")

load("results_robust4_continue.RData")
ssx4_continue = results$ssx
ssx4_continue = data.frame(ssx4_continue)
ssx4_continue$method = "pilot + summ"
colnames(ssx4_continue) = c("s4", "method")



ssx_results = rbind(ssx,ssx_pilot,ssx4,ssx4_continue)

p4 = ggplot(data=ssx_results) + 
  geom_density(aes(x=s4, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position = c(0.3, 0.85)) +
  theme(legend.title=element_blank()) +
  geom_point(aes(x=obs_summ[4], y=0), colour="black", size = 4) +
  xlim(c(1.72,1.77)) + theme(axis.title=element_text(size=16)) + labs(y= "density")
 


x11()
grid.arrange(p1, p2, p3, p4, nrow = 1)







