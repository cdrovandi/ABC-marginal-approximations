smc_abc_generic_continue <- function(theta,ssx,N,a,acc_rate_stop,prior_sim,prior_eval,distance_fun,distance_fun_pilot,extra_args){
  
 # generic implementation of the SMC ABC replenishment algorithm of Drovandi and Pettitt (2011) Biometrics.
 # this version of the code is initialised with an ABC sample rather than a prior sample.
  
  param = prior_sim(extra_args)
  num_params = length(param)
  
  N_keep = ceiling(a*N)
  N_drop = N-N_keep
  
  trial_mcmc_iters = 5
  
  dist = rep(0,N)
 
  
  for (i in 1:N){
    dist[i] = sum((extra_args$obs_summ - ssx[i,])^2)
  }
  
  
  acc_rate = 1
  
  while (acc_rate > acc_rate_stop){
    
    sort_dist = sort(dist, index.return = TRUE)
    dist = sort_dist$x
    theta = theta[sort_dist$ix,]
    if (extra_args$return_summ == TRUE){
      ssx = as.matrix(ssx[sort_dist$ix,])
    }
    
    dist_next = dist[N_keep] # adaptively determine next distance
    sprintf("the next ABC tolerance is %s", dist_next)
    
    # resampling
    ind = sample(1:N_keep, N_drop, replace = TRUE)
    
    dist[(N_keep+1):N] = dist[ind]
    theta[(N_keep+1):N,] = theta[ind,]
    if (extra_args$return_summ == TRUE){
      ssx[(N_keep+1):N,] = ssx[ind,]
    }
    
    cov_rw = cov(theta)
    accept_count = rep(0,N_drop)
    
    cat("number of trial mcmc iterations is ",trial_mcmc_iters, "\n")
    for (i in (N_keep+1):N){
      for (r in 1:trial_mcmc_iters){
        theta_prop = mvrnorm(n=1, mu = theta[i,], Sigma = cov_rw)
        prior_curr = prior_eval(theta[i,], extra_args)
        prior_prop = prior_eval(theta_prop, extra_args)
        if (runif(1) > prior_prop/prior_curr){  # early rejection
          next
        }
        ret = distance_fun_pilot(theta_prop, extra_args)
        dist_prop = ret$dist
        if (dist_prop <= dist_next){
          theta[i,] = theta_prop
          dist[i] = dist_prop
          if (extra_args$return_summ == TRUE){
            ssx[i,] = ret$sim_summ  
          }
          accept_count[i-N_keep] = accept_count[i-N_keep] + 1
        }
      }
    }
    
    acc_rate = sum(accept_count)/(N_drop*trial_mcmc_iters)
    mcmc_iters = ceiling(log(0.01)/log(1-acc_rate))
    remaining_mcmc_iters = mcmc_iters - trial_mcmc_iters
    
    cat("number of remaining mcmc iterations is ",remaining_mcmc_iters, "\n")
    for (i in (N_keep+1):N){
      for (r in 1:remaining_mcmc_iters){
        theta_prop = mvrnorm(n=1, mu = theta[i,], Sigma = cov_rw)
        prior_curr = prior_eval(theta[i,], extra_args)
        prior_prop = prior_eval(theta_prop, extra_args)
        if (runif(1) > prior_prop/prior_curr){  # early rejection
          next
        }
        ret = distance_fun_pilot(theta_prop, extra_args)
        dist_prop = ret$dist
        if (dist_prop <= dist_next){
          theta[i,] = theta_prop
          dist[i] = dist_prop
          if (extra_args$return_summ == TRUE){
            ssx[i,] = ret$sim_summ  
          }
          accept_count[i-N_keep] = accept_count[i-N_keep] + 1
        }
      }
    }
    
    acc_rate = sum(accept_count)/(N_drop*mcmc_iters)
    cat("MCMC acceptance rate was ", acc_rate, "\n")
    cat("The number of unique particles is ", length(unique(theta[,1])), "\n")
    
    trial_mcmc_iters = ceiling(0.5*mcmc_iters)
    
  }
  
  results = list()
  
  results$theta = theta
  results$dist = dist
  if (extra_args$return_summ == TRUE){
    results$ssx = ssx 
  }
  
  return(results)
  
  
  
}


