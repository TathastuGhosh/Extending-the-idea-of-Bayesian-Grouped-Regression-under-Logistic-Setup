
# ================================================================================
# Generating beta for each setting
# ================================================================================

get_beta_and_gs = function(setting){
  if(setting == 1){
    beta = numeric(50)
    pos1 = sample(1:10, 1)
    pos2 = 10 + sample(1:10, 1)
    pos3 = 20 + sample(1:10, 1)
    pos4 = 30 + sample(1:10, 1)
    pos5 = 40 + sample(1:10, 1)
    beta[c(pos1,pos2,pos3,pos4,pos5)] = c(0.5,1,1.5,2,2)
    return(list(beta = beta, group_sizes = rep(10,5)))
  }
  if(setting == 2){
    beta = numeric(50)
    beta[1:10] = 0.5
    return(list(beta = beta, group_sizes = rep(10,5)))
  }
  if(setting == 3){
    beta = c(
      {x=rep(0,30); x[sample(30,27)] = sample(c(3,5,8),27,TRUE); x},
      {x=rep(0,10); x[sample(10,8)] = sample(c(2,3,4,5),8,TRUE); x},
      {x=rep(0,20); x[sample(20,18)] = 7.5; x},
      {x=rep(0,5);  x[sample(5,4)]  = sample(c(9,5,8,7),4); x},
      {x=rep(0,15); x[sample(15,14)] = 1.5; x},
      {x=rep(0,20); x[sample(20,18)] = 0.5; x}
    )
    return(list(beta = beta, group_sizes = c(30,10,20,5,15,20)))
  }
  if(setting == 4){
    beta = c(
      {x = rep(0,25); x[sample(25,22)] = 0.8; x},
      {x = rep(0,10);x},
      {x = rep(0,10); x[sample(10,3)] = 2.5; x},
      {x = rep(0,10); x[sample(10,8)] = 1.5; x},
      {x = rep(0, 5); x},
      {x = rep(0,15); x[sample(15,4)] = sample(c(1,2,3,5), 4); x},
      {x = rep(0,5); x[sample(5,1)] = 2; x}
    )
    return(list(beta = beta,group_sizes = c(25,10,10,10,5,15,5)))
  }
}


# =====================================================================
# Generating data for logistic grouped setup
# =====================================================================

generate_data = function(n, rho_w= 0.8,rho_b = 0.2, setting){
  
  if(!(setting %in% 1:4)) 
    stop("Setting should be 1,2,3 or 4.")
  
  output = get_beta_and_gs(setting)
  beta = output$beta
  group_sizes = output$group_sizes
  
  
  # Correlation Matrix
  groups = rep(seq_along(group_sizes), group_sizes)
  Sigma = outer(groups, groups,
                FUN = function(x, y) ifelse(x == y, rho_w, rho_b))
  
  diag(Sigma) = 1
  if(any(eigen(Sigma)$values <= 0))
    stop("Sigma is not positive definite.")
  
  # Generate X
  X = MASS::mvrnorm(n = n, mu = rep(0, length(beta)), Sigma = Sigma)
  
  # Generate probabilities
  eta = X %*% beta
  pi = plogis(eta)
  
  # Generating Bernoulli response
  y = rbinom(n, size = 1, prob = pi)
  
  # Final data frame
  data = data.frame(y, X)
  colnames(data) = c("y", paste0("x",1:ncol(X)))
  
  return(list(data = data, beta = beta, 
              Sigma = Sigma, group_sizes = group_sizes))
}


generate_data(n = 10, setting = 1)








