


# Package names
packages <- c("MASS", "MCMCpack", "tidyverse", "mvtnorm", "invgamma", "ggplot2", "gridExtra", "devtools")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

install_github("jrlewi/brlm")
library(brlm)

# load libraries
#library(MASS)
#library(brlm)
#library(MCMCpack)
#library(tidyverse)
#library(mvtnorm)
#library(invgamma)
#library(ggplot2)
#library(gridExtra)



source("auxiliary_functions.R")
source("smc_abc_generic.R")
source("smc_abc_generic_continue.R")


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





############ SMC ABC with all summaries ###############

 
extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ = compute_statistics(y, extra_args)
extra_args$obs_summ = obs_summ

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance,extra_args)

save(file="results_summ.RData", results)



############ SMC ABC with all summaries (pilot) ###############


extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ = compute_statistics(y, extra_args)
extra_args$obs_summ = obs_summ

results = smc_abc_generic(1000,0.5,0.25,prior_sim_trans,prior_eval,distance,extra_args)

save(file="results_summ_pilot.RData", results)


############ SMC ABC with first summary ###############


extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ = compute_statistics(y, extra_args)
obs_summ = obs_summ[1]
extra_args$obs_summ = obs_summ

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance1,extra_args)

save(file="results_summ1.RData", results)


############ SMC ABC with second summary ###############


extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ = compute_statistics(y, extra_args)
obs_summ = obs_summ[2]
extra_args$obs_summ = obs_summ

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance2,extra_args)

save(file="results_summ2.RData", results)


############ SMC ABC with third summary ###############


extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ = compute_statistics(y, extra_args)
obs_summ = obs_summ[3]
extra_args$obs_summ = obs_summ

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance3,extra_args)

save(file="results_summ3.RData", results)



############ SMC ABC with fourth summary ###############


extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ = compute_statistics(y, extra_args)
obs_summ = obs_summ[4]
extra_args$obs_summ = obs_summ

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance4,extra_args)

save(file="results_summ4.RData", results)



############ SMC ABC with first summary using pilot run ###############


load("results_summ_pilot.RData")
extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ_pilot = compute_statistics(y, extra_args)
obs_summ = compute_statistics1(y, extra_args)
dist_pilot = max(results$dist)

extra_args$obs_summ_pilot = obs_summ_pilot
extra_args$obs_summ = obs_summ
extra_args$dist_pilot = dist_pilot

results = smc_abc_generic_continue(results$theta,as.matrix(results$ssx[,1]),1000,0.5,0.01,prior_sim_trans,prior_eval,distance1,distance_pilot1,extra_args)

save(results, file = "results_summ1_continue.RData")



############ SMC ABC with second summary using pilot run ###############


load("results_summ_pilot.RData")
extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ_pilot = compute_statistics(y, extra_args)
obs_summ = compute_statistics2(y, extra_args)
dist_pilot = max(results$dist)

extra_args$obs_summ_pilot = obs_summ_pilot
extra_args$obs_summ = obs_summ
extra_args$dist_pilot = dist_pilot

results = smc_abc_generic_continue(results$theta,as.matrix(results$ssx[,2]),1000,0.5,0.01,prior_sim_trans,prior_eval,distance2,distance_pilot2,extra_args)

save(results, file = "results_summ2_continue.RData")


############ SMC ABC with third summary using pilot run ###############


load("results_summ_pilot.RData")
extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ_pilot = compute_statistics(y, extra_args)
obs_summ = compute_statistics3(y, extra_args)
dist_pilot = max(results$dist)

extra_args$obs_summ_pilot = obs_summ_pilot
extra_args$obs_summ = obs_summ
extra_args$dist_pilot = dist_pilot

results = smc_abc_generic_continue(results$theta,as.matrix(results$ssx[,3]),1000,0.5,0.01,prior_sim_trans,prior_eval,distance3,distance_pilot3,extra_args)

save(results, file = "results_summ3_continue.RData")





############ SMC ABC with fourth summary using pilot run ###############


load("results_summ_pilot.RData")
extra_args = list(X = X, psi = psi.huber, scaleEst = 'Huber', maxit = 100, N=N, p=p, return_summ = TRUE)
obs_summ_pilot = compute_statistics(y, extra_args)
obs_summ = compute_statistics4(y, extra_args)
dist_pilot = max(results$dist)

extra_args$obs_summ_pilot = obs_summ_pilot
extra_args$obs_summ = obs_summ
extra_args$dist_pilot = dist_pilot

results = smc_abc_generic_continue(results$theta,as.matrix(results$ssx[,4]),1000,0.5,0.01,prior_sim_trans,prior_eval,distance4,distance_pilot4,extra_args)

save(results, file = "results_summ4_continue.RData")




