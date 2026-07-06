

# MH for a_g^2
# -----------------------------

mh_update_a_g_sq_gigg = function(a_g_2_current, delta_g_sq, d_g, sigma_prop){
  
  log_a_g_2_new = rnorm(1,log(a_g_2_current), sigma_prop)
  a_g_2_new = exp(log_a_g_2_new)
  
  
  log_alpha = -(lgamma(sqrt(a_g_2_new)) - lgamma(sqrt(a_g_2_current))) +
    (sqrt(a_g_2_new)-sqrt(a_g_2_current)) * log(delta_g_sq) -
    0.5*(log(a_g_2_new)-log(a_g_2_current)) -
    (1/d_g)*((1/a_g_2_new) - (1/a_g_2_current))
  log_alpha = min(0, log_alpha)
  
  if(log(runif(1)) <= log_alpha){
    return(a_g_2_new)
  }else{
    return(a_g_2_current)
  }
}


# MH for b_g^2
# ----------------------------------

mh_update_b_g_sq_gigg = function(b_g_2_current, lambda_vals, e_g, p_g, sigma_prop){
  
  S_t = sum(log(lambda_vals))
  log_b_g_2_new = rnorm(1, log(b_g_2_current), sigma_prop)
  b_g_2_new = exp(log_b_g_2_new)
  
  log_alpha = p_g*(lgamma(sqrt(b_g_2_current))-lgamma(sqrt(b_g_2_new))) -
    (sqrt(b_g_2_new)-sqrt(b_g_2_current))*S_t -
    0.5*(log(b_g_2_new)-log(b_g_2_current)) -
    (1/e_g)*((1/b_g_2_new)-(1/b_g_2_current))
  log_alpha = min(0,log_alpha)
  
  if(log(runif(1)) <= log_alpha){
    return(b_g_2_new)
  }else{
    return(b_g_2_current)
  }
}





# GIGG Gibbs
# ----------------------------------------

gigg_gibbs = function(data, nsamples = 1000, burnin = 1000,
                      init, group_sizes, print_every = 100,
                      sigma_prop_a = 0.2, sigma_prop_b = 0.2){
  require(MASS); require(BayesLogit); require(MCMCpack); require(GIGrvg)
  
  n = nrow(data)
  p = ncol(data) - 1
  X = as.matrix(data[,-1])
  y = as.vector(data[,1])
  G = length(group_sizes)
  
  expected_len = n + p + p + G + 1 + 1 + 4*G
  if(length(init) != expected_len) stop("init length does not match number of parameters.")
  if(burnin < 0) stop("burnin must be non-negative.")
  if(burnin >= nsamples) stop("burnin must be smaller than nsamples.")
  if(sum(group_sizes) != p) stop("sum(group_sizes) must equal number of covariates.")
  
  # creating the sample dataframe
  colname = c(paste0("w_",1:n), paste0("beta_",1:p),
              unlist(lapply(1:G, function(g) paste0("lam_sq_", g, "_", 1:group_sizes[g]))),
              paste0("delta_sq_",1:G),
              "tau_sq","v",
              paste0("a_sq_",1:G), paste0("d_",1:G), 
              paste0("b_sq_",1:G), paste0("e_",1:G)
              )
  sample = matrix(NA, nrow = nsamples, ncol = length(colname))
  colnames(sample) = colname
  
  # extract initial values
  idx = 0
  current_w = init[(idx+1):(idx+n)]; idx = idx + n
  current_beta = init[(idx+1):(idx+p)]; idx = idx + p
  current_lambda_sq = init[(idx+1):(idx+p)]; idx = idx + p
  current_delta_sq = init[(idx+1):(idx+G)]; idx = idx + G
  current_tau_sq = init[idx+1]; idx = idx + 1
  current_v = init[idx+1]; idx = idx + 1
  current_a_sq = init[(idx+1):(idx+G)]; idx = idx + G
  current_d = init[(idx+1):(idx+G)]; idx = idx + G
  current_b_sq = init[(idx+1):(idx+G)]; idx = idx + G
  current_e = init[(idx+1):(idx+G)]; idx = idx + G
  
  idx_groups = split(1:p, rep(1:G, times = group_sizes))
  k = y - 0.5
  
  # initial row
  sample[1,] = c(current_w, current_beta, current_lambda_sq, 
                 current_delta_sq, current_tau_sq, current_v, current_a_sq,
                 current_d, current_b_sq, current_e)
  start_time = Sys.time()
  cat("\n===============================================================\n")
  cat("                Starting GIGG Gibbs Sampler\n")
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
    if(iter %% print_every == 0){
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      avg_time <- elapsed/(iter-1)
      eta_remaining <- avg_time * (nsamples - iter)
      pct <- round(100 * iter/nsamples, 1)
      current_time <- format(Sys.time(), "%H:%M:%S")
      cat(sprintf("Time: %s => Samples: [%d/%d] | Completed: %.1f%% | elapsed: %.1fs | ETA: %.1fs\n",
                  current_time, iter, nsamples, pct, elapsed, eta_remaining))
    }
    
    # ---------- 1. Sample w ----------------
    eta = as.vector(X %*% current_beta)
    current_w = BayesLogit::rpg(num = n, h = 1, z = eta)
    
    # ---------- 2. Sample beta --------------
    D_vec = current_tau_sq * current_delta_sq[rep(1:G, times = group_sizes)] * current_lambda_sq
    D_inv = diag(1/D_vec)
    XtOmegaX = crossprod(X, X * current_w)
    A = XtOmegaX + D_inv
    R = chol(A)
    Sigma_beta = chol2inv(R)
    mu_beta = Sigma_beta%*% crossprod(X,k)
    current_beta = MASS::mvrnorm(n = 1, mu = as.vector(mu_beta),
                                 Sigma = Sigma_beta)
    # ----------- 3. Sample lambda sq -----------
    b_g_vec = sqrt(current_b_sq[rep(1:G, times = group_sizes)])
    delta_g_vec = rep(current_delta_sq, times = group_sizes)
    shape_lambda = b_g_vec + 0.5
    scale_lambda = 1 + (current_beta^2)/(2*current_tau_sq*delta_g_vec)
    current_lambda_sq = rinvgamma(p, shape = shape_lambda, scale = scale_lambda)
    
    # ------------ 4. Sample delta sq -----------
    new_delta_sq = numeric(G)
    a_g_vec = sqrt(current_a_sq)
    for(g in 1:G){
      ind = idx_groups[[g]]
      p_g = group_sizes[g]
      sum_term = sum((current_beta[ind]^2)/(current_tau_sq * current_lambda_sq[ind]))
      lambda_gig = a_g_vec[g] - p_g/2
      chi_gig = sum_term
      psi_gig = 2.0
      new_delta_sq[g] = rgig(1, lambda = lambda_gig, chi = chi_gig, psi = psi_gig)
    }
    current_delta_sq = new_delta_sq
    
    # ------------ 5. Sample tau_sq and v --------
    total_sum = sum((current_beta^2) / 
                      (current_lambda_sq * rep(current_delta_sq, times = group_sizes)))
    shape_tau = (p + 1)/2
    scale_tau = 0.5*total_sum + 1/current_v
    current_tau_sq = rinvgamma(1, shape = shape_tau, scale = scale_tau)
    current_v = rinvgamma(1, shape = 1, scale = 1 + 1/current_tau_sq)
    
    # ------------ 6. Auxiliary variables  -------
    current_d = rinvgamma(G, shape = 1, scale = 1+1/current_a_sq)
    current_e = rinvgamma(G, shape = 1, scale = 1+1/current_b_sq)
    
    # ------------ 7. Metropolis - Hatings update for a_sq and b_sq -----
    # a_g^2 (per group)
    for(g in 1:G){
      current_a_sq[g] = mh_update_a_g_sq_gigg(current_a_sq[g],current_delta_sq[g],
                                              current_d[g], sigma_prop_a)
    }
    # b_g^2 (per group)
    for(g in 1:G){
      ind = idx_groups[[g]]
      lambda_vals = current_lambda_sq[ind]
      p_g = group_sizes[g]
      current_b_sq[g] = mh_update_b_g_sq_gigg(current_b_sq[g], lambda_vals,
                                              current_e[g], p_g, sigma_prop_b)
    }
    
    # -------------- 8. Store Samples ------------------
    new_row = c(current_w, current_beta, current_lambda_sq,
                current_delta_sq, current_tau_sq, current_v,
                current_a_sq, current_d, current_b_sq, current_e)
    sample[iter,] = new_row
  }
  
  # Discard burn-in
  sample = as.data.frame(sample)
  post_idx = seq(burnin+1, nsamples, by = 1)
  posterior_sample = sample[post_idx,]
  rownames(posterior_sample) = seq_len(nrow(posterior_sample))
  
  elapsed = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
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






# ------------------------------------------
# Checking
# ------------------------------------------

dat = generate_data(n = 500, setting = 1)
group_sizes = dat$group_sizes
n = nrow(dat$data)
p = ncol(dat$data) - 1
G = length(group_sizes)

create_gigg_init = function(n,p,G, group_sizes){
  c(rep(0.5,n),            # w
    rep(0,p),              # beta
    rep(1,p),              # lambda_sq
    rep(1,G),              # delta_sq
    1,1,                   # tau_sq and v
    rep(1,G),              # a_sq
    rep(1,G),              # d
    rep(1,G),              # b_sq
    rep(1,G)               # e
  )
}

init = create_gigg_init(n,p,G, group_sizes)
post_gigg = gigg_gibbs(data = dat$data, nsamples = 50000, burnin = 10000,
                       init = init, group_sizes = group_sizes,
                       sigma_prop_a = 0.2, sigma_prop_b = 0.2)
























































































































