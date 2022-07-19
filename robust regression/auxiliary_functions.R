



compute_statistics <- function(y, extra_args){
  
  X = extra_args$X
  psi = extra_args$psi
  scaleEst = extra_args$scaleEst
  maxit = extra_args$maxit
  
  robust <- MASS::rlm(X, y, psi = psi, scale.est = scaleEst, maxit = maxit)
  c(robust$coef, log(robust$s))
}


compute_statistics1 <- function(y, extra_args){
  
  X = extra_args$X
  psi = extra_args$psi
  scaleEst = extra_args$scaleEst
  maxit = extra_args$maxit
  
  robust <- MASS::rlm(X, y, psi = psi, scale.est = scaleEst, maxit = maxit)
  robust$coef[1]
}

compute_statistics2 <- function(y, extra_args){
  
  X = extra_args$X
  psi = extra_args$psi
  scaleEst = extra_args$scaleEst
  maxit = extra_args$maxit
  
  robust <- MASS::rlm(X, y, psi = psi, scale.est = scaleEst, maxit = maxit)
  robust$coef[2]
}


compute_statistics3 <- function(y, extra_args){
  
  X = extra_args$X
  psi = extra_args$psi
  scaleEst = extra_args$scaleEst
  maxit = extra_args$maxit
  
  robust <- MASS::rlm(X, y, psi = psi, scale.est = scaleEst, maxit = maxit)
  robust$coef[3]
}

compute_statistics4 <- function(y, extra_args){
  
  X = extra_args$X
  psi = extra_args$psi
  scaleEst = extra_args$scaleEst
  maxit = extra_args$maxit
  
  robust <- MASS::rlm(X, y, psi = psi, scale.est = scaleEst, maxit = maxit)
  log(robust$s)
}



# function to simulate the prior (here we transform parameters to the real line)
prior_sim_trans <- function(extra_args){

  return(c(rnorm(3,0,10), log(rexp(1))))
  
}


# function to evaluate prior density (for transformed parameters)
prior_eval <- function(theta, extra_args){
  
  return(prod(c(dnorm(theta[1:3],0,10), dexp(exp(theta[4]),1)*exp(theta[4]))))

}




# distance function calculation using all summaries
distance <- function(theta, extra_args){
  
  N = extra_args$N
  X = extra_args$X
  p = extra_args$p
  ys <- rnorm(N, X%*%as.matrix(theta[1:p]), exp(theta[p+1]))
  
  
  sim_summ = compute_statistics(ys, extra_args)
  dist <- sum((extra_args$obs_summ - sim_summ)^2)
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}


# distance function calculation using first summary
distance1 <- function(theta, extra_args){
  
  N = extra_args$N
  X = extra_args$X
  p = extra_args$p
  ys <- rnorm(N, X%*%as.matrix(theta[1:p]), exp(theta[p+1]))
  
  
  sim_summ = compute_statistics(ys, extra_args)
  sim_summ = sim_summ[1]
  dist <- sum((extra_args$obs_summ - sim_summ)^2)
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}

# distance function calculation using second summary
distance2 <- function(theta, extra_args){
  
  N = extra_args$N
  X = extra_args$X
  p = extra_args$p
  ys <- rnorm(N, X%*%as.matrix(theta[1:p]), exp(theta[p+1]))
  
  
  sim_summ = compute_statistics(ys, extra_args)
  sim_summ = sim_summ[2]
  dist <- sum((extra_args$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}


# distance function calculation using third summary
distance3 <- function(theta, extra_args){
  
  N = extra_args$N
  X = extra_args$X
  p = extra_args$p
  ys <- rnorm(N, X%*%as.matrix(theta[1:p]), exp(theta[p+1]))
  
  
  sim_summ = compute_statistics(ys, extra_args)
  sim_summ = sim_summ[3]
  dist <- sum((extra_args$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}


# distance function calculation using fourth summary
distance4 <- function(theta, extra_args){
  
  N = extra_args$N
  X = extra_args$X
  p = extra_args$p
  ys <- rnorm(N, X%*%as.matrix(theta[1:p]), exp(theta[p+1]))
  
  
  sim_summ = compute_statistics(ys, extra_args)
  sim_summ = sim_summ[4]
  dist <- sum((extra_args$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}



# distance function calculation using first summary including pilot distance
distance_pilot1 <- function(theta, inp){
  
  
  N = extra_args$N
  X = extra_args$X
  p = extra_args$p
  ys <- rnorm(N, X%*%as.matrix(theta[1:p]), exp(theta[p+1]))
  
  sim_summ = compute_statistics(ys, extra_args)
  
  dist_prop <- sum((inp$obs_summ_pilot - sim_summ)^2)
  
  if (dist_prop > inp$dist_pilot){ # return very large distance (it will get rejected)
    ret = list()
    ret$dist = 10000
    return(ret)
  }
  
  sim_summ = compute_statistics1(ys, extra_args)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}




# distance function calculation using second summary including pilot distance
distance_pilot2 <- function(theta, inp){
  
  
  N = extra_args$N
  X = extra_args$X
  p = extra_args$p
  ys <- rnorm(N, X%*%as.matrix(theta[1:p]), exp(theta[p+1]))
  
  sim_summ = compute_statistics(ys, extra_args)
  
  dist_prop <- sum((inp$obs_summ_pilot - sim_summ)^2)
  
  if (dist_prop > inp$dist_pilot){ # return very large distance (it will get rejected)
    ret = list()
    ret$dist = 10000
    return(ret)
  }
  
  sim_summ = compute_statistics2(ys, extra_args)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}




# distance function calculation using third summary including pilot distance
distance_pilot3 <- function(theta, inp){
  
  
  N = extra_args$N
  X = extra_args$X
  p = extra_args$p
  ys <- rnorm(N, X%*%as.matrix(theta[1:p]), exp(theta[p+1]))
  
  sim_summ = compute_statistics(ys, extra_args)
  
  dist_prop <- sum((inp$obs_summ_pilot - sim_summ)^2)
  
  if (dist_prop > inp$dist_pilot){ # return very large distance (it will get rejected)
    ret = list()
    ret$dist = 10000
    return(ret)
  }
  
  sim_summ = compute_statistics3(ys, extra_args)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}






# distance function calculation using fourth summary including pilot distance
distance_pilot4 <- function(theta, inp){
  
  
  N = extra_args$N
  X = extra_args$X
  p = extra_args$p
  ys <- rnorm(N, X%*%as.matrix(theta[1:p]), exp(theta[p+1]))
  
  sim_summ = compute_statistics(ys, extra_args)
  
  dist_prop <- sum((inp$obs_summ_pilot - sim_summ)^2)
  
  if (dist_prop > inp$dist_pilot){ # return very large distance (it will get rejected)
    ret = list()
    ret$dist = 10000
    return(ret)
  }
  
  sim_summ = compute_statistics4(ys, extra_args)
  
  dist <- sum((inp$obs_summ - sim_summ)^2)
  
  
  ret = list()
  
  ret$dist = dist
  if (extra_args$return_summ == TRUE){
    ret$sim_summ = sim_summ
  }
  
  return(ret)
}




