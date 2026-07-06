
# MH for a_g^2
# ----------------------------------

mh_update_a_g_sq_grasp = function(a_g_2_current, t_vals, d_g, p_g, sigma_prop) {
  S_t = sum(log(t_vals))
  log_a_g_2_new = rnorm(1, log(a_g_2_current), sigma_prop)
  a_g_2_new = exp(log_a_g_2_new)
  
  log_alpha = -p_g * (lgamma(sqrt(a_g_2_new)) - lgamma(sqrt(a_g_2_current))) -
    (sqrt(a_g_2_new) - sqrt(a_g_2_current)) * S_t -
    0.5 * (log(a_g_2_new) - log(a_g_2_current)) -
    (1 / d_g) * ((1 / a_g_2_new) - (1 / a_g_2_current))
  log_alpha = min(0, log_alpha)
  
  if (log(runif(1)) <= log_alpha) {
    return(a_g_2_new)
  } else {
    return(a_g_2_current)
  }
}


# MH for b_g^2
# # ----------------------------------

mh_update_b_g_sq_grasp = function(b_g_2_current, lambda_2_vals, t_vals, e_g, p_g, sigma_prop){
  S_t = sum(log(t_vals))
  S_lambda = sum(log(lambda_2_vals))
  
  log_b_g_2_new = rnorm(1, log(b_g_2_current), sigma_prop)
  b_g_2_new = exp(log_b_g_2_new)
  
  log_alpha = -p_g*(lgamma(sqrt(b_g_2_new)) - lgamma(sqrt(b_g_2_current))) -
    (sqrt(b_g_2_new)-sqrt(b_g_2_current))*S_t -
    (sqrt(b_g_2_new)-sqrt(b_g_2_current))*S_lambda -
    0.5*(log(b_g_2_new)-log(b_g_2_current)) -
    (1/e_g)*((1/b_g_2_new) - (1/b_g_2_current))
  log_alpha = min(0, log_alpha)
  
  if (log(runif(1)) <= log_alpha) {
    return(b_g_2_new)
  } else {
    return(b_g_2_current)
  }
}



# MH for a^2
# ----------------------------------

mh_update_a_sq_grasp = function(a_2_current, c_vals, d, G, sigma_prop){
  S_c = sum(log(c_vals))
  
  log_a_2_new = rnorm(1, log(a_2_current), sigma_prop)
  a_2_new = exp(log_a_2_new)
  
  log_alpha = -G*(lgamma(sqrt(a_2_new)) - lgamma(sqrt(a_2_current))) - 
    (sqrt(a_2_new)-sqrt(a_2_current))*S_c -
    0.5*(log(a_2_new) - log(a_2_current)) -
    (1/d) * ((1/a_2_new) - (1/a_2_current))
  log_alpha = min(0, log_alpha)
  
  if(log(runif(1))<= log_alpha){
    return(a_2_new)
  }else{
    return(a_2_current)
  }
}

# MH for b^2
# ----------------------------------

mh_update_b_sq_grasp = function(b_2_current, delta_2_vals, c_vals, e, G, sigma_prop){
  
  S_delta = sum(log(delta_2_vals))
  S_c = sum(log(c_vals))
  log_b_2_new = rnorm(1, log(b_2_current), sigma_prop)
  b_2_new = exp(log_b_2_new)
  
  
  log_alpha = -G * (lgamma(sqrt(b_2_new))-lgamma(sqrt(b_2_current))) -
    (sqrt(b_2_new) - sqrt(b_2_current))*S_c - (sqrt(b_2_new) - sqrt(b_2_current))*S_delta -
    0.5*(log(b_2_new) - log(b_2_current)) -
    (1/e) * ((1/b_2_new) - (1/b_2_current))
  log_alpha = min(0, log_alpha)
  
  if(log(runif(1))<= log_alpha){
    return(b_2_new)
  }else{
    return(b_2_current)
  }
}





# GRASP Gibbs
# ----------------------------------

grasp_gibbs = function(data, nsamples = 10000, burnin = 1000,
                       init, group_sizes, print_every = 100,
                       sigma_prop_a = 0.2, sigma_prop_b = 0.2,
                       sigma_prop_a_global = 0.2, sigma_prop_b_global = 0.2){
  
  require(MASS)
  require(BayesLogit)
  require(MCMCpack)
  
  n = nrow(data)
  p = ncol(data) - 1
  X = as.matrix(data[,-1])
  y = as.vector(data[,1])
  G = length(group_sizes)
  
  expected_len = n + 3*p + 2*G + 2 + 4*G + 4
  if(length(init) != expected_len) stop("init length does not match number of parameters.")
  if(burnin < 0) stop("burnin must be non-negative.")
  if(burnin >= nsamples) stop("burnin must be smaller than nsamples.")
  if(sum(group_sizes) != p) stop("sum(group_sizes) must equal number of covariates.")
  
  # Creating the sample dataframe
  colname = c(paste0("w_",1:n),paste0("beta_",1:p),
              unlist(Map(function(i,n) paste0("lam_sq_", i, 1:n), 1:G, group_sizes)),
              unlist(Map(function(i,n) paste0("t_", i, 1:n), 1:G, group_sizes)),
              paste0("delta_sq_",1:G),
              paste0("c_",1:G),
              "tau_sq","v",
              paste0("a_sq_",1:G),
              paste0("d_",1:G),
              paste0("b_sq_",1:G),
              paste0("e_",1:G),
              "a_sq","d","b_sq","e"
              )
  
  sample= matrix(NA, nrow = nsamples, ncol = length(colname))
  colnames(sample) = colname
  
  # Tweaking initials
  idx = 0
  current_w = init[(idx+1):(idx+n)]; idx = idx + n
  current_beta = init[(idx+1):(idx+p)]; idx = idx + p
  current_lambda_sq = init[(idx+1):(idx+p)]; idx = idx + p
  current_t = init[(idx+1):(idx+p)]; idx = idx + p
  current_delta_sq = init[(idx+1):(idx+G)]; idx = idx + G
  current_c = init[(idx+1):(idx+G)]; idx = idx + G
  current_tau_sq = init[idx+1]; idx = idx + 1
  current_v = init[idx+1]; idx = idx + 1
  current_a_g_sq = init[(idx+1):(idx+G)]; idx = idx + G
  current_d_g = init[(idx+1):(idx+G)]; idx = idx + G
  current_b_g_sq = init[(idx+1):(idx+G)]; idx = idx + G
  current_e_g = init[(idx+1):(idx+G)]; idx = idx + G
  current_a_sq_global = init[idx+1]; idx = idx + 1
  current_d_global = init[idx+1]; idx = idx + 1
  current_b_sq_global = init[idx+1]; idx = idx + 1
  current_e_global = init[idx+1]
  
  # Few predefined things
  idx_groups = split(1:p, rep(1:G, times = group_sizes))
  k = y - 0.5
  
  # initial row 
  sample[1,] = c(current_w, current_beta, current_lambda_sq, current_t,
                 current_delta_sq, current_c, current_tau_sq, current_v,
                 current_a_g_sq, current_d_g, current_b_g_sq, current_e_g,
                 current_a_sq_global, current_d_global, current_b_sq_global, 
                 current_e_global)
  
  # Time related (Initial)
  start_time = Sys.time()
  cat("\n===============================================================\n")
  cat("                Starting GRASP Gibbs Sampler\n")
  cat("===============================================================\n\n")
  cat("----------------------------------------\n")
  cat("Summary of experimental setting\n")
  cat("----------------------------------------\n")
  cat(sprintf("Observations    : %d\n", n))
  cat(sprintf("Covariates      : %d\n", p))
  cat(sprintf("Number of groups: %d\n", G))
  cat(sprintf("Group sizes     : %s\n", paste(group_sizes, collapse = ", ")))
  cat("\n----------------------------------------\n\n")
  
  for(iter in 2:nsamples){
    
    # Progress Printing
    if(iter %% print_every == 0){
      elapsed = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      avg_time = elapsed/(iter-1)
      eta_remaining = avg_time * (nsamples - iter)
      pct = round(100 * iter/nsamples,1)
      current_time = format(Sys.time(),"%H:%M:%S")
      cat(sprintf("Time: %s => Samples: [%d/%d] | Completed: %.1f%% | elapsed: %.1fs | ETA: %.1fs\n",
                  current_time, iter, nsamples, pct, elapsed, eta_remaining))
    }
    
    # --------------- 1. Sample w ---------------------
    eta = as.vector(X %*% current_beta)
    current_w = BayesLogit::rpg(num = n, h = 1, z = eta)
    
    # --------------- 2. Sample beta ------------------
    D_vec = current_tau_sq * current_lambda_sq * rep(current_delta_sq, times = group_sizes)
    D_inv = diag(1/D_vec)
    XtOmegaX = crossprod(X, X * current_w)
    P = XtOmegaX + D_inv
    R = chol(P)
    mu_beta = backsolve(R, forwardsolve(t(R), crossprod(X,k)))
    z = rnorm(p)
    current_beta = mu_beta + backsolve(R, z)
    
    # ---------------- 3. Sample lambda ---------------
    b_g_vec = sqrt(current_b_g_sq[rep(1:G, times = group_sizes)])
    delta_g_vec = rep(current_delta_sq, times = group_sizes)
    scale_lambda = (current_beta^2) / (2*current_tau_sq*delta_g_vec) + 1/current_t
    current_lambda_sq = rinvgamma(p, shape = b_g_vec + 1/2, scale = scale_lambda)
    
    # --------------- 4. Sample t --------------------
    a_g_vec = sqrt(current_a_g_sq[rep(1:G, times = group_sizes)])
    shape_t = a_g_vec + b_g_vec
    scale_t = 1+(1/(current_lambda_sq))
    current_t = rinvgamma(p, shape = shape_t, scale = scale_t)
    
    # 5. Sample delta_sq
    b_global = sqrt(current_b_sq_global)
    new_delta_sq = numeric(G)
    for (g in 1:G) {
      ind = idx_groups[[g]]
      shape_delta = b_global + group_sizes[g] / 2
      sum_term = sum(current_beta[ind]^2 / current_lambda_sq[ind])
      scale_delta = sum_term / (2 * current_tau_sq) + 1 / current_c[g]
      new_delta_sq[g] = rinvgamma(1, shape = shape_delta, scale = scale_delta)
    }
    current_delta_sq = new_delta_sq
    
    # 6. Sample c
    a_global = sqrt(current_a_sq_global)
    shape_c = a_global + b_global
    scale_c = 1 + (1/current_delta_sq)
    current_c = rinvgamma(G, shape = shape_c, scale = scale_c)
    
    # 7. Sample tau_Sq and v
    total_sum = sum(current_beta^2 / (current_lambda_sq * rep(current_delta_sq, times = group_sizes)))
    shape_tau = (p + 1) / 2
    scale_tau = 0.5 * total_sum + 1 / current_v
    current_tau_sq = rinvgamma(1, shape = shape_tau, scale = scale_tau)
    current_v = rinvgamma(1, shape = 1, scale = 1 + 1 / current_tau_sq)
    
    
    # 8. Sample auxiliary variables for hyperpriors (conjugate)
    current_d_g <- rinvgamma(G, shape = 1, scale = 1 + 1 / current_a_g_sq)
    current_e_g <- rinvgamma(G, shape = 1, scale = 1 + 1 / current_b_g_sq)
    current_d_global <- rinvgamma(1, shape = 1, scale = 1 + 1 / current_a_sq_global)
    current_e_global <- rinvgamma(1, shape = 1, scale = 1 + 1 / current_b_sq_global)
    
    # 9. Metropolis-Hastings updates for squared hyperparameters
    # a_g_sq (per group)
    for (g in 1:G) {
      ind <- idx_groups[[g]]
      t_vals <- current_t[ind]
      current_a_g_sq[g] <- mh_update_a_g_sq_grasp(current_a_g_sq[g], t_vals,
                                            current_d_g[g], group_sizes[g], sigma_prop_a)
    }
    # b_g_sq (per group)
    for (g in 1:G) {
      ind <- idx_groups[[g]]
      lambda_vals <- current_lambda_sq[ind]
      t_vals <- current_t[ind]
      current_b_g_sq[g] <- mh_update_b_g_sq_grasp(current_b_g_sq[g], lambda_vals, t_vals,
                                            current_e_g[g], group_sizes[g], sigma_prop_b)
    }
    # a_sq_global
    current_a_sq_global <- mh_update_a_sq_grasp(current_a_sq_global, current_c,
                                          current_d_global, G, sigma_prop_a_global)
    # b_sq_global
    current_b_sq_global <- mh_update_b_sq_grasp(current_b_sq_global, current_delta_sq, current_c,
                                          current_e_global, G, sigma_prop_b_global)
    
    # 10. Store samples
    new_row <- c(current_w, current_beta, current_lambda_sq, current_t,
                 current_delta_sq, current_c, current_tau_sq, current_v,
                 current_a_g_sq, current_d_g, current_b_g_sq, current_e_g,
                 current_a_sq_global, current_d_global, current_b_sq_global, current_e_global)
    sample[iter, ] <- new_row
  }
  
  # Discard burn‑in and return
  sample <- as.data.frame(sample)
  post_idx <- seq(burnin + 1, nsamples, by = 1)
  posterior_sample <- sample[post_idx, ]
  rownames(posterior_sample) <- seq_len(nrow(posterior_sample))
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat("\n----------------------------------------\n")
  cat("Run Summary\n")
  cat("----------------------------------------\n")
  cat(sprintf("Total samples : %d\n", nsamples))
  cat(sprintf("Burn-in       : %d\n", burnin))
  cat(sprintf("Posterior kept: %d\n", nrow(posterior_sample)))
  elapsed_min <- floor(elapsed / 60)
  elapsed_sec <- round(elapsed %% 60)
  cat(sprintf("Elapsed time  : %.1f seconds (%dm %02ds)\n", elapsed, elapsed_min, elapsed_sec))
  cat("----------------------------------------\n\n")
  
  return(posterior_sample)
  
}



# -----------------------------------
# Checking
# -----------------------------------

dat = generate_data(n = 500, setting = 1)
group_sizes = dat$group_sizes
n = nrow(dat$data)
p =ncol(dat$data) - 1
G = length(group_sizes)

create_grasp_init = function(n, p, G, group_sizes) {
  c(rep(0.5, n),                # w
    rep(0, p),                  # beta
    rep(1, p),                  # lambda_sq
    rep(1, p),                  # t
    rep(1, G),                  # delta_sq
    rep(1, G),                  # c
    1,                          # tau_sq
    1,                          # v
    rep(1, G),                  # a_g_sq
    rep(1, G),                  # d_g
    rep(1, G),                  # b_g_sq
    rep(1, G),                  # e_g
    1,                          # a_sq_global
    1,                          # d_global
    1,                          # b_sq_global
    1)                          # e_global
}

init = create_grasp_init(n,p,G,group_sizes)
post = grasp_gibbs(data = dat$data, nsamples = 10000, burnin = 0,
                   init = init, group_sizes = group_sizes,
                   sigma_prop_a = 0.2, sigma_prop_b = 0.2,
                   sigma_prop_a_global = 0.2, sigma_prop_b_global = 0.2)
View(post)
















































































































