library(ggplot2)
library(gridExtra)
library(grid)

library(MASS)
library(brlm)
library(MCMCpack)
library(tidyverse)
library(mvtnorm)
library(invgamma)
library(ggplot2)
library(gridExtra)

source("auxiliary_functions.R")

State_keep = 27

trend_prior <- sqrt_count_2010 ~ sqrt_count_2008 - 1 + 
  Associate_Count +
  Office_Employee_Count

trend <- sqrt_count_2012 ~ sqrt_count_2010 - 1 + 
  Associate_Count +
  Office_Employee_Count

prior_data <- read_rds('prior_data.rds')


prior_data <- prior_data %>% 
  mutate(sqrt_count_2008 = sqrt(Count_2008), 
         sqrt_count_2010 = sqrt(Count_2010))  %>% 
  filter(State == State_keep, Type == 1, Count_2010 > 0)


# pooled regression analysis ----
prior_fit <- MASS::rlm(trend_prior, 
                       scale.est = 'Huber', 
                       data =  prior_data, maxit = 100)

analysis_data <- read_rds('analysis_data.rds')

analysis_data <- analysis_data %>% 
  mutate(sqrt_count_2010 = sqrt(Count_2010), 
         sqrt_count_2012 = sqrt(Count_2012)) %>% 
  filter(State == State_keep)

N <- nrow(analysis_data)
p <- length(coef(prior_fit))

# prior hyperparameters
beta_0 <- coef(prior_fit)
sigma_beta_0 <- vcov(prior_fit)
var_scalar <- floor(nrow(prior_data))
var_beta_0 <- var_scalar*sigma_beta_0
sigma2_hat <- prior_fit$s^2

a_0 <- 5
b_0 <- sigma2_hat*(a_0 - 1)   


X = cbind(analysis_data$sqrt_count_2010, analysis_data$Associate_Count, analysis_data$Office_Employee_Count)

y = analysis_data$sqrt_count_2012

extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ = compute_statistics(y, extra_args)


################## parameter 1 ##############


load("results_summ.RData")
theta = results$theta
theta[,4] = exp(theta[,4])^2
theta = data.frame(theta)
theta$method = "all summs"
colnames(theta) = c("b0", "b1", "b2", "sigma^2", "method")

load("results_summ_pilot.RData")
theta_pilot = results$theta
theta_pilot[,4] = exp(theta_pilot[,4])^2
theta_pilot = data.frame(theta_pilot)
theta_pilot$method = "pilot (all summs)"
colnames(theta_pilot) = c("b0", "b1", "b2", "sigma^2", "method")

load("results_summ1.RData")
theta1 = results$theta
theta1[,4] = exp(theta1[,4])^2
theta1 = data.frame(theta1)
theta1$method = "summ 1"
colnames(theta1) = c("b0", "b1", "b2", "sigma^2", "method")

load("results_summ1_continue.RData")
theta1_continue = results$theta
theta1_continue[,4] = exp(theta1_continue[,4])^2
theta1_continue = data.frame(theta1_continue)
theta1_continue$method = "pilot + summ 1"
colnames(theta1_continue) = c("b0", "b1", "b2", "sigma^2", "method")



theta_results = rbind(theta,theta_pilot,theta1,theta1_continue)

p1 = ggplot(data=theta_results) + 
  geom_density(aes(x=b0, group = method, colour = method, linetype = method), size = 1.5) +
  #theme(legend.position = c(0.2, 0.85)) +
  #theme(legend.title=element_blank()) +
#+
  theme(legend.position="none") +
 # geom_point(aes(x=3, y=0), colour="gold4", size = 4) +
  xlim(c(0.75,1.25)) + theme(axis.title=element_text(size=16))





################## parameter 2 ##############


load("results_summ.RData")
theta = results$theta
theta[,4] = exp(theta[,4])^2
theta = data.frame(theta)
theta$method = "all summs"
colnames(theta) = c("b0", "b1", "b2", "sigma^2", "method")

load("results_summ_pilot.RData")
theta_pilot = results$theta
theta_pilot[,4] = exp(theta_pilot[,4])^2
theta_pilot = data.frame(theta_pilot)
theta_pilot$method = "pilot (all summs)"
colnames(theta_pilot) = c("b0", "b1", "b2", "sigma^2", "method")

load("results_summ2.RData")
theta2 = results$theta
theta2[,4] = exp(theta2[,4])^2
theta2 = data.frame(theta2)
theta2$method = "summ 2"
colnames(theta2) = c("b0", "b1", "b2", "sigma^2", "method")

load("results_summ2_continue.RData")
theta2_continue = results$theta
theta2_continue[,4] = exp(theta2_continue[,4])^2
theta2_continue = data.frame(theta2_continue)
theta2_continue$method = "pilot + summ 2"
colnames(theta2_continue) = c("b0", "b1", "b2", "sigma^2", "method")



theta_results = rbind(theta,theta_pilot,theta2,theta2_continue)

p2 = ggplot(data=theta_results) + 
  geom_density(aes(x=b1, group = method, colour = method, linetype = method), size = 1.5) +
  #theme(legend.position = c(0.2, 0.85)) +
  #theme(legend.title=element_blank()) +
#+
  theme(legend.position="none") +
# geom_point(aes(x=3, y=0), colour="gold4", size = 4) +
  xlim(c(-0.2,0.2)) + theme(axis.title=element_text(size=16))






################## parameter 3 ##############


load("results_summ.RData")
theta = results$theta
theta[,4] = exp(theta[,4])^2
theta = data.frame(theta)
theta$method = "all summs"
colnames(theta) = c("b0", "b1", "b2", "sigma^2", "method")

load("results_summ_pilot.RData")
theta_pilot = results$theta
theta_pilot[,4] = exp(theta_pilot[,4])^2
theta_pilot = data.frame(theta_pilot)
theta_pilot$method = "pilot (all summs)"
colnames(theta_pilot) = c("b0", "b1", "b2", "sigma^2", "method")

load("results_summ3.RData")
theta3 = results$theta
theta3[,4] = exp(theta3[,4])^2
theta3 = data.frame(theta3)
theta3$method = "summ 3"
colnames(theta3) = c("b0", "b1", "b2", "sigma^2", "method")

load("results_summ3_continue.RData")
theta3_continue = results$theta
theta3_continue[,4] = exp(theta3_continue[,4])^2
theta3_continue = data.frame(theta3_continue)
theta3_continue$method = "pilot + summ 3"
colnames(theta3_continue) = c("b0", "b1", "b2", "sigma^2", "method")



theta_results = rbind(theta,theta_pilot,theta3,theta3_continue)

p3 = ggplot(data=theta_results) + 
  geom_density(aes(x=b2, group = method, colour = method, linetype = method), size = 1.5) +
  #theme(legend.position = c(0.2, 0.85)) +
  #theme(legend.title=element_blank()) +
  theme(legend.position="none") +
#+
# geom_point(aes(x=3, y=0), colour="gold4", size = 4) +
  xlim(c(-0.25,0.25)) + theme(axis.title=element_text(size=16))
 







################## parameter 4 ##############


load("results_summ.RData")
theta = results$theta
theta[,4] = exp(theta[,4])^2
theta = data.frame(theta)
theta$method = "all summs"
colnames(theta) = c("b0", "b1", "b2", "sigma2", "method")

load("results_summ_pilot.RData")
theta_pilot = results$theta
theta_pilot[,4] = exp(theta_pilot[,4])^2
theta_pilot = data.frame(theta_pilot)
theta_pilot$method = "pilot (all summs)"
colnames(theta_pilot) = c("b0", "b1", "b2", "sigma2", "method")

load("results_summ4.RData")
theta4 = results$theta
theta4[,4] = exp(theta4[,4])^2
theta4 = data.frame(theta4)
theta4$method = "summ"
colnames(theta4) = c("b0", "b1", "b2", "sigma2", "method")

load("results_summ4_continue.RData")
theta4_continue = results$theta
theta4_continue[,4] = exp(theta4_continue[,4])^2
theta4_continue = data.frame(theta4_continue)
theta4_continue$method = "pilot + summ"
colnames(theta4_continue) = c("b0", "b1", "b2", "sigma2", "method")



theta_results = rbind(theta,theta_pilot,theta4,theta4_continue)

p4 = ggplot(data=theta_results) + 
  geom_density(aes(x=sigma2, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position = c(0.7, 0.85)) +
  theme(legend.title=element_blank()) + theme(axis.title=element_text(size=16))
#+
# geom_point(aes(x=3, y=0), colour="gold4", size = 4) +
#  xlim(c(2.5,3.5))




x11()
grid.arrange(p1, p2, p3, p4, nrow = 1)








################## s1 ##############


load("results_summ.RData")
ssx = results$ssx[,1]
ssx = data.frame(ssx)
ssx$method = "all summs"
colnames(ssx) = c("s1", "method")

load("results_summ_pilot.RData")
ssx_pilot = results$ssx[,1]
ssx_pilot = data.frame(ssx_pilot)
ssx_pilot$method = "pilot (all summs)"
colnames(ssx_pilot) = c("s1", "method")

load("results_summ1.RData")
ssx1 = results$ssx
ssx1 = data.frame(ssx1)
ssx1$method = "summ 1"
colnames(ssx1) = c("s1", "method")

load("results_summ1_continue.RData")
ssx1_continue = results$ssx
ssx1_continue = data.frame(ssx1_continue)
ssx1_continue$method = "pilot + summ 1"
colnames(ssx1_continue) = c("s1", "method")



ssx_results = rbind(ssx,ssx_pilot,ssx1,ssx1_continue)

p1 = ggplot(data=ssx_results) + 
  geom_density(aes(x=s1, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position="none") +
  geom_point(aes(x=obs_summ[1], y=0), colour="black", size = 4) +
  xlim(c(0.95,0.99)) + theme(axis.title=element_text(size=16)) + labs(y= "density")



################## s2 ##############


load("results_summ.RData")
ssx = results$ssx[,2]
ssx = data.frame(ssx)
ssx$method = "all summs"
colnames(ssx) = c("s2", "method")

load("results_summ_pilot.RData")
ssx_pilot = results$ssx[,2]
ssx_pilot = data.frame(ssx_pilot)
ssx_pilot$method = "pilot (all summs)"
colnames(ssx_pilot) = c("s2", "method")

load("results_summ2.RData")
ssx2 = results$ssx
ssx2 = data.frame(ssx2)
ssx2$method = "summ 2"
colnames(ssx2) = c("s2", "method")

load("results_summ2_continue.RData")
ssx2_continue = results$ssx
ssx2_continue = data.frame(ssx2_continue)
ssx2_continue$method = "pilot + summ 2"
colnames(ssx2_continue) = c("s2", "method")



ssx_results = rbind(ssx,ssx_pilot,ssx2,ssx2_continue)

p2 = ggplot(data=ssx_results) + 
  geom_density(aes(x=s2, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position="none") +
  geom_point(aes(x=obs_summ[2], y=0), colour="black", size = 4) +
  xlim(c(-0.01,0.03)) + theme(axis.title=element_text(size=16)) + labs(y= "density")





################## s3 ##############


load("results_summ.RData")
ssx = results$ssx[,3]
ssx = data.frame(ssx)
ssx$method = "all summs"
colnames(ssx) = c("s3", "method")

load("results_summ_pilot.RData")
ssx_pilot = results$ssx[,3]
ssx_pilot = data.frame(ssx_pilot)
ssx_pilot$method = "pilot (all summs)"
colnames(ssx_pilot) = c("s3", "method")

load("results_summ3.RData")
ssx3 = results$ssx
ssx3 = data.frame(ssx3)
ssx3$method = "summ 3"
colnames(ssx3) = c("s3", "method")

load("results_summ3_continue.RData")
ssx3_continue = results$ssx
ssx3_continue = data.frame(ssx3_continue)
ssx3_continue$method = "pilot + summ 3"
colnames(ssx3_continue) = c("s3", "method")



ssx_results = rbind(ssx,ssx_pilot,ssx3,ssx3_continue)

p3 = ggplot(data=ssx_results) + 
  geom_density(aes(x=s3, group = method, colour = method, linetype = method), size = 1.5) +
  theme(legend.position="none") +
  geom_point(aes(x=obs_summ[3], y=0), colour="black", size = 4) +
  xlim(c(0.01,0.06)) + theme(axis.title=element_text(size=16)) + labs(y= "density")






################## s4 ##############


load("results_summ.RData")
ssx = results$ssx[,4]
ssx = data.frame(ssx)
ssx$method = "all summs"
colnames(ssx) = c("s4", "method")

load("results_summ_pilot.RData")
ssx_pilot = results$ssx[,4]
ssx_pilot = data.frame(ssx_pilot)
ssx_pilot$method = "pilot (all summs)"
colnames(ssx_pilot) = c("s4", "method")

load("results_summ4.RData")
ssx4 = results$ssx
ssx4 = data.frame(ssx4)
ssx4$method = "summ"
colnames(ssx4) = c("s4", "method")

load("results_summ4_continue.RData")
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
  xlim(c(-2.55,-2.4)) + theme(axis.title=element_text(size=16)) + labs(y= "density")



x11()
grid.arrange(p1, p2, p3, p4, nrow = 1)


