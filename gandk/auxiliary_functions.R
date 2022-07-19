
# auxiliary functions required by SMC ABC


# function to simulate the prior (here we transform parameters to the real line)
prior_sim_trans <- function(extra_args){
  
  lower = extra_args$lower
  upper = extra_args$upper
  
  theta = runif(length(lower), lower, upper)
  return(log((theta - lower)/(upper - theta)))
  
  
}


# function to simulate the prior (here we transform parameters to the real line)
prior_sim <- function(extra_args){
  
  lower = extra_args$lower
  upper = extra_args$upper
  
  theta = runif(length(lower), lower, upper)
  return(theta)
  
  
}


# function to evaluate prior density (for transformed parameters)
prior_eval <- function(theta, extra_args){
  
  lower = extra_args$lower
  upper = extra_args$upper
  
  prior_value <- prod((upper - lower)*exp(theta)/((1+exp(theta))^2))
  
  return(prior_value)
}



# summary statistic function for robust summaries
summary_robust <- function(x){
  octiles = quantile(x,prob=seq(from=1/8,to=7/8,length=7))
  momb = octiles[6]-octiles[2]
  c(octiles[4], momb, (octiles[6]+octiles[2]-2*octiles[4])/momb, (octiles[7]-octiles[5]+octiles[3] - octiles[1])/momb)
}

# summary statistic function for first robust summary
summary_robust1 <- function(x){
  octiles = quantile(x,prob=seq(from=1/8,to=7/8,length=7))
  momb = octiles[6]-octiles[2]
  c(octiles[4])
}

# summary statistic function for second robust summary
summary_robust2 <- function(x){
  octiles = quantile(x,prob=seq(from=1/8,to=7/8,length=7))
  momb = octiles[6]-octiles[2]
  c(momb)
}


# summary statistic function for third robust summary
summary_robust3 <- function(x){
  octiles = quantile(x,prob=seq(from=1/8,to=7/8,length=7))
  momb = octiles[6]-octiles[2]
  c((octiles[6]+octiles[2]-2*octiles[4])/momb)
}


# summary statistic function for fourth robust summary
summary_robust4 <- function(x){
  octiles = quantile(x,prob=seq(from=1/8,to=7/8,length=7))
  momb = octiles[6]-octiles[2]
  c((octiles[7]-octiles[5]+octiles[3] - octiles[1])/momb)
}


# distance function calculation using 4 robust summaries
distance_robust <- function(theta, extra_args){
  
  theta = (extra_args$lower*exp(-theta) + extra_args$upper)/(exp(-theta) + 1)
  
  x <- rgk(extra_args$nobs,theta[1],theta[2],theta[3],theta[4])
  sim_summ = summary_robust(x)
  
  dist <- sum((extra_args$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}


# distance function calculation using first robust summary
distance_robust1 <- function(theta, inp){
  
  theta = (inp$lower*exp(-theta) + inp$upper)/(exp(-theta) + 1)
  
  x <- rgk(inp$nobs,theta[1],theta[2],theta[3],theta[4])
  sim_summ = summary_robust1(x)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (inp$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}



# distance function calculation using second robust summary
distance_robust2 <- function(theta, inp){
  
  theta = (inp$lower*exp(-theta) + inp$upper)/(exp(-theta) + 1)
  
  x <- rgk(inp$nobs,theta[1],theta[2],theta[3],theta[4])
  sim_summ = summary_robust2(x)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  ret = list()
  
  ret$dist = dist
  if (inp$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}


# distance function calculation using third robust summary
distance_robust3 <- function(theta, inp){
  
  theta = (inp$lower*exp(-theta) + inp$upper)/(exp(-theta) + 1)
  
  x <- rgk(inp$nobs,theta[1],theta[2],theta[3],theta[4])
  sim_summ = summary_robust3(x)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (inp$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}


# distance function calculation using fourth robust summary
distance_robust4 <- function(theta, inp){
  
  theta = (inp$lower*exp(-theta) + inp$upper)/(exp(-theta) + 1)
  
  x <- rgk(inp$nobs,theta[1],theta[2],theta[3],theta[4])
  sim_summ = summary_robust4(x)
  
  output <- sum((inp$obs_summ - sim_summ)^2)
  
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (inp$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}




# distance function calculation using first robust summary including pilot distance
distance_pilot_robust1 <- function(theta, inp){
  
  theta = (inp$lower*exp(-theta) + inp$upper)/(exp(-theta) + 1)
  
  x <- rgk(inp$nobs,theta[1],theta[2],theta[3],theta[4])
  sim_summ = summary_robust(x)
  
  dist_prop <- sum((inp$obs_summ_pilot - sim_summ)^2)
  
  if (dist_prop > inp$dist_pilot){ # return very large distance
    ret = list()
    ret$dist = 10000
    return(ret)
  }
  
  sim_summ = summary_robust1(x)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (inp$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}


# distance function calculation using second robust summary including pilot distance
distance_pilot_robust2 <- function(theta, inp){
  
  theta = (inp$lower*exp(-theta) + inp$upper)/(exp(-theta) + 1)
  
  x <- rgk(inp$nobs,theta[1],theta[2],theta[3],theta[4])
  sim_summ = summary_robust(x)
  
  dist_prop <- sum((inp$obs_summ_pilot - sim_summ)^2)
  
  if (dist_prop > inp$dist_pilot){ # return very large distance
    ret = list()
    ret$dist = 10000
    return(ret)
  }
  
  sim_summ = summary_robust2(x)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (inp$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}


# distance function calculation using third robust summary including pilot distance
distance_pilot_robust3 <- function(theta, inp){
  
  theta = (inp$lower*exp(-theta) + inp$upper)/(exp(-theta) + 1)
  
  x <- rgk(inp$nobs,theta[1],theta[2],theta[3],theta[4])
  sim_summ = summary_robust(x)
  
  dist_prop <- sum((inp$obs_summ_pilot - sim_summ)^2)
  
  if (dist_prop > inp$dist_pilot){ # return very large distance
    ret = list()
    ret$dist = 10000
    return(ret)
  }
  
  sim_summ = summary_robust3(x)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (inp$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}


# distance function calculation using fourth robust summary including pilot distance
distance_pilot_robust4 <- function(theta, inp){
  
  theta = (inp$lower*exp(-theta) + inp$upper)/(exp(-theta) + 1)
  
  x <- rgk(inp$nobs,theta[1],theta[2],theta[3],theta[4])
  sim_summ = summary_robust(x)
  
  dist_prop <- sum((inp$obs_summ_pilot - sim_summ)^2)
  
  if (dist_prop > inp$dist_pilot){ # return very large distance
    ret = list()
    ret$dist = 10000
    return(ret)
  }
  
  sim_summ = summary_robust4(x)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (inp$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}





