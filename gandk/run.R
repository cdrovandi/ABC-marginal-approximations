


# libraries
library(parallel)
library(gk)
library(mvtnorm)
library(MASS)

 
source("smc_abc_generic.R")
source("smc_abc_generic_continue.R")
source("auxiliary_functions.R")

# generate observed data
#y = rgk(10000,3,1,2,0.5)
#save(y, file = "data_gandk.RData")

# load the "observed" data
load(file = "data_gandk.RData")


############ SMC ABC with robust summaries ###############

obs_summ = summary_robust(y)
extra_args = list(obs_summ = obs_summ, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance_robust,extra_args)
save(results, file = "results_robust.RData")

############ SMC ABC with robust summaries pilot run ###############

obs_summ = summary_robust(y)
extra_args = list(obs_summ = obs_summ, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic(1000,0.5,0.2,prior_sim_trans,prior_eval,distance_robust,extra_args)
save(results, file = "results_robust_pilot.RData")

############ SMC ABC with first robust summary ###############

obs_summ = summary_robust1(y)
extra_args = list(obs_summ = obs_summ, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance_robust1,extra_args)
save(results, file = "results_robust1.RData")



############ SMC ABC with second robust summary ###############

obs_summ = summary_robust2(y)
extra_args = list(obs_summ = obs_summ, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance_robust2,extra_args)
save(results, file = "results_robust2.RData")


############ SMC ABC with third robust summary ###############

obs_summ = summary_robust3(y)
extra_args = list(obs_summ = obs_summ, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance_robust3,extra_args)
save(results, file = "results_robust3.RData")


############ SMC ABC with fourth robust summary ###############

obs_summ = summary_robust4(y)
extra_args = list(obs_summ = obs_summ, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic(1000,0.5,0.01,prior_sim_trans,prior_eval,distance_robust4,extra_args)
save(results, file = "results_robust4.RData")


############ SMC ABC with first robust summary using pilot run ###############


load("results_robust_pilot.RData")
obs_summ_pilot = summary_robust(y)
obs_summ = summary_robust1(y)
dist_pilot = max(results$dist)

extra_args = list(obs_summ = obs_summ, obs_summ_pilot = obs_summ_pilot, dist_pilot = dist_pilot, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic_continue(results$theta,as.matrix(results$ssx[,1]),1000,0.5,0.01,prior_sim_trans,prior_eval,distance_robust1,distance_pilot_robust1,extra_args)

save(results, file = "results_robust1_continue.RData")


############ SMC ABC with second robust summary using pilot run ###############


load("results_robust_pilot.RData")
obs_summ_pilot = summary_robust(y)
obs_summ = summary_robust2(y)
dist_pilot = max(results$dist)

extra_args = list(obs_summ = obs_summ, obs_summ_pilot = obs_summ_pilot, dist_pilot = dist_pilot, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic_continue(results$theta,as.matrix(results$ssx[,2]),1000,0.5,0.01,prior_sim_trans,prior_eval,distance_robust2,distance_pilot_robust2,extra_args)

save(results, file = "results_robust2_continue.RData")


############ SMC ABC with third robust summary using pilot run ###############


load("results_robust_pilot.RData")
obs_summ_pilot = summary_robust(y)
obs_summ = summary_robust3(y)
dist_pilot = max(results$dist)

extra_args = list(obs_summ = obs_summ, obs_summ_pilot = obs_summ_pilot, dist_pilot = dist_pilot, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic_continue(results$theta,as.matrix(results$ssx[,3]),1000,0.5,0.01,prior_sim_trans,prior_eval,distance_robust3,distance_pilot_robust3,extra_args)

save(results, file = "results_robust3_continue.RData")


############ SMC ABC with fourth robust summary using pilot run ###############


load("results_robust_pilot.RData")
obs_summ_pilot = summary_robust(y)
obs_summ = summary_robust4(y)
dist_pilot = max(results$dist)

extra_args = list(obs_summ = obs_summ, obs_summ_pilot = obs_summ_pilot, dist_pilot = dist_pilot, nobs = 10000, lower = c(0, 0, 0, 0), upper = c(10, 10, 10, 10), return_summ = TRUE)

results = smc_abc_generic_continue(results$theta,as.matrix(results$ssx[,4]),1000,0.5,0.01,prior_sim_trans,prior_eval,distance_robust4,distance_pilot_robust4,extra_args)

save(results, file = "results_robust4_continue.RData")


