# ================================================================
# GRASP Gibbs Sampler (with optional return of acceptance rates)
# ================================================================


grasp_gibbs_for_tuning <- function(data, nsamples = 10000, burnin = 1000,
                                   init, group_sizes, print_every = 100,
                                   sigma_prop_a = 0.2, sigma_prop_b = 0.2,
                                   sigma_prop_a_global = 0.2, sigma_prop_b_global = 0.2,
                                   return_accept = FALSE) {
  
  require(MASS); require(BayesLogit); require(MCMCpack)
  
  n <- nrow(data)
  p <- ncol(data) - 1
  X <- as.matrix(data[,-1])
  y <- as.vector(data[,1])
  G <- length(group_sizes)
  
  expected_len <- n + 3*p + 2*G + 2 + 4*G + 4
  if(length(init) != expected_len) stop("init length mismatch")
  if(burnin < 0) stop("burnin non‑negative")
  if(burnin >= nsamples) stop("burnin < nsamples")
  if(sum(group_sizes) != p) stop("group_sizes sum != p")
  
  # Column names (same as in grasp_gibbs)
  colname <- c(paste0("w_",1:n), paste0("beta_",1:p),
               unlist(Map(function(i,n) paste0("lam_sq_", i, 1:n), 1:G, group_sizes)),
               unlist(Map(function(i,n) paste0("t_", i, 1:n), 1:G, group_sizes)),
               paste0("delta_sq_",1:G), paste0("c_",1:G),
               "tau_sq","v",
               paste0("a_sq_",1:G), paste0("d_",1:G),
               paste0("b_sq_",1:G), paste0("e_",1:G),
               "a_sq","d","b_sq","e")
  
  sample <- matrix(NA, nrow = nsamples, ncol = length(colname))
  colnames(sample) <- colname
  
  # Extract initial values
  idx <- 0
  current_w <- init[(idx+1):(idx+n)]; idx <- idx + n
  current_beta <- init[(idx+1):(idx+p)]; idx <- idx + p
  current_lambda_sq <- init[(idx+1):(idx+p)]; idx <- idx + p
  current_t <- init[(idx+1):(idx+p)]; idx <- idx + p
  current_delta_sq <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_c <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_tau_sq <- init[idx+1]; idx <- idx + 1
  current_v <- init[idx+1]; idx <- idx + 1
  current_a_g_sq <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_d_g <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_b_g_sq <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_e_g <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_a_sq_global <- init[idx+1]; idx <- idx + 1
  current_d_global <- init[idx+1]; idx <- idx + 1
  current_b_sq_global <- init[idx+1]; idx <- idx + 1
  current_e_global <- init[idx+1]
  
  idx_groups <- split(1:p, rep(1:G, times = group_sizes))
  k <- y - 0.5
  
  # Acceptance counters (only used if return_accept = TRUE)
  if (return_accept) {
    accept_a <- rep(0, G)
    accept_b <- rep(0, G)
    accept_a_global <- 0
    accept_b_global <- 0
    total_post_mh <- 0
  }
  
  # Initial row
  sample[1,] <- c(current_w, current_beta, current_lambda_sq, current_t,
                  current_delta_sq, current_c, current_tau_sq, current_v,
                  current_a_g_sq, current_d_g, current_b_g_sq, current_e_g,
                  current_a_sq_global, current_d_global, current_b_sq_global, current_e_global)
  
  start_time <- Sys.time()
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
  
  for(iter in 2:nsamples) {
    if(iter %% print_every == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
      avg_time <- elapsed/(iter-1)
      eta <- avg_time * (nsamples - iter)
      pct <- round(100 * iter/nsamples, 1)
      cat(sprintf("Iter %d/%d (%.1f%%) | elapsed %.1fs | ETA %.1fs\n",
                  iter, nsamples, pct, elapsed, eta))
    }
    
    # ---------- 1. Sample w ----------
    eta <- as.vector(X %*% current_beta)
    current_w <- BayesLogit::rpg(num = n, h = 1, z = eta)
    
    # ---------- 2. Sample beta ----------
    D_vec <- current_tau_sq * current_lambda_sq * rep(current_delta_sq, times = group_sizes)
    D_inv <- diag(1/D_vec)
    XtOmegaX <- crossprod(X, X * current_w)
    P <- XtOmegaX + D_inv
    R <- chol(P)
    mu_beta <- backsolve(R, forwardsolve(t(R), crossprod(X, k)))
    current_beta <- mu_beta + backsolve(R, rnorm(p))
    
    # ---------- 3. Sample lambda_sq ----------
    b_g_vec <- sqrt(current_b_g_sq[rep(1:G, times = group_sizes)])
    delta_g_vec <- rep(current_delta_sq, times = group_sizes)
    scale_lambda <- (current_beta^2) / (2*current_tau_sq*delta_g_vec) + 1/current_t
    current_lambda_sq <- rinvgamma(p, shape = b_g_vec + 0.5, scale = scale_lambda)
    
    # ---------- 4. Sample t ----------
    a_g_vec <- sqrt(current_a_g_sq[rep(1:G, times = group_sizes)])
    shape_t <- a_g_vec + b_g_vec
    scale_t <- 1 + 1/current_lambda_sq
    current_t <- rinvgamma(p, shape = shape_t, scale = scale_t)
    
    # ---------- 5. Sample delta_sq ----------
    b_global <- sqrt(current_b_sq_global)
    new_delta_sq <- numeric(G)
    for (g in 1:G) {
      ind <- idx_groups[[g]]
      shape_delta <- b_global + group_sizes[g]/2
      sum_term <- sum(current_beta[ind]^2 / current_lambda_sq[ind])
      scale_delta <- sum_term/(2*current_tau_sq) + 1/current_c[g]
      new_delta_sq[g] <- rinvgamma(1, shape = shape_delta, scale = scale_delta)
    }
    current_delta_sq <- new_delta_sq
    
    # ---------- 6. Sample c ----------
    a_global <- sqrt(current_a_sq_global)
    shape_c <- a_global + b_global
    scale_c <- 1 + 1/current_delta_sq
    current_c <- rinvgamma(G, shape = shape_c, scale = scale_c)
    
    # ---------- 7. Sample tau and v ----------
    total_sum <- sum(current_beta^2 / (current_lambda_sq * rep(current_delta_sq, times = group_sizes)))
    current_tau_sq <- rinvgamma(1, shape = (p+1)/2, scale = 0.5*total_sum + 1/current_v)
    current_v <- rinvgamma(1, shape = 1, scale = 1 + 1/current_tau_sq)
    
    # ---------- 8. Auxiliaries for hyperpriors ----------
    current_d_g <- rinvgamma(G, shape = 1, scale = 1 + 1/current_a_g_sq)
    current_e_g <- rinvgamma(G, shape = 1, scale = 1 + 1/current_b_g_sq)
    current_d_global <- rinvgamma(1, shape = 1, scale = 1 + 1/current_a_sq_global)
    current_e_global <- rinvgamma(1, shape = 1, scale = 1 + 1/current_b_sq_global)
    
    # ---------- 9. MH updates with acceptance counting ----------
    # a_g_sq
    for (g in 1:G) {
      ind <- idx_groups[[g]]
      t_vals <- current_t[ind]
      old <- current_a_g_sq[g]
      new <- mh_update_a_g_sq_grasp(old, t_vals, current_d_g[g], group_sizes[g], sigma_prop_a)
      if (return_accept && iter > burnin && new != old) accept_a[g] <- accept_a[g] + 1
      current_a_g_sq[g] <- new
    }
    # b_g_sq
    for (g in 1:G) {
      ind <- idx_groups[[g]]
      lambda_vals <- current_lambda_sq[ind]
      t_vals <- current_t[ind]
      old <- current_b_g_sq[g]
      new <- mh_update_b_g_sq_grasp(old, lambda_vals, t_vals, current_e_g[g], group_sizes[g], sigma_prop_b)
      if (return_accept && iter > burnin && new != old) accept_b[g] <- accept_b[g] + 1
      current_b_g_sq[g] <- new
    }
    # a_sq_global
    old <- current_a_sq_global
    new <- mh_update_a_sq_grasp(old, current_c, current_d_global, G, sigma_prop_a_global)
    if (return_accept && iter > burnin && new != old) accept_a_global <- accept_a_global + 1
    current_a_sq_global <- new
    # b_sq_global
    old <- current_b_sq_global
    new <- mh_update_b_sq_grasp(old, current_delta_sq, current_c, current_e_global, G, sigma_prop_b_global)
    if (return_accept && iter > burnin && new != old) accept_b_global <- accept_b_global + 1
    current_b_sq_global <- new
    
    if (return_accept && iter > burnin) total_post_mh <- total_post_mh + 1
    
    # ---------- 10. Store samples ----------
    sample[iter,] <- c(current_w, current_beta, current_lambda_sq, current_t,
                       current_delta_sq, current_c, current_tau_sq, current_v,
                       current_a_g_sq, current_d_g, current_b_g_sq, current_e_g,
                       current_a_sq_global, current_d_global, current_b_sq_global, current_e_global)
  }
  
  # Discard burn‑in
  sample <- as.data.frame(sample)
  post_idx <- seq(burnin+1, nsamples)
  posterior_sample <- sample[post_idx,]
  rownames(posterior_sample) <- seq_len(nrow(posterior_sample))
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
  cat("\n----------------------------------------\n")
  cat("Run Summary\n")
  cat("----------------------------------------\n")
  cat(sprintf("Total samples: %d\n", nsamples))
  cat(sprintf("Burn-in: %d\n", burnin))
  cat(sprintf("Posterior kept: %d\n", nrow(posterior_sample)))
  cat(sprintf("Elapsed time: %.1f seconds\n", elapsed))
  cat("----------------------------------------\n\n")
  
  if (return_accept && total_post_mh > 0) {
    acceptance <- list(
      a_g = accept_a / total_post_mh,
      b_g = accept_b / total_post_mh,
      a_global = accept_a_global / total_post_mh,
      b_global = accept_b_global / total_post_mh
    )
    return(list(posterior = posterior_sample, acceptance = acceptance))
  } else if (return_accept) {
    warning("No post‑burn‑in MH steps; acceptance rates NA")
    return(list(posterior = posterior_sample, acceptance = list(a_g=NA, b_g=NA, a_global=NA, b_global=NA)))
  } else {
    return(posterior_sample)
  }
}
# ============================================================
# Tune proposal sigmas for GRASP using a grid search
# ============================================================
tune_sigma_grasp <- function(data, init, group_sizes,
                             sigma_candidates = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0),
                             target_rate = 0.3,
                             pilot_nsamples = 800,
                             pilot_burnin = 300,
                             print_every = 200) {
  
  n <- nrow(data)
  p <- ncol(data) - 1
  G <- length(group_sizes)
  
  # Data frame to store results
  results <- data.frame(sigma = sigma_candidates,
                        acc_a = NA, acc_b = NA,
                        acc_a_global = NA, acc_b_global = NA)
  
  cat("\n=============================================================\n")
  cat("Tuning proposal sigma for GRASP\n")
  cat("=============================================================\n")
  cat(sprintf("Target acceptance rate: %.2f\n", target_rate))
  cat(sprintf("Pilot: %d iterations, burn‑in %d\n", pilot_nsamples, pilot_burnin))
  cat("Testing sigmas:", paste(sigma_candidates, collapse=", "), "\n")
  cat("=============================================================\n\n")
  
  for (s in sigma_candidates) {
    cat(sprintf("\nTesting sigma = %.3f ...\n", s))
    
    pilot <- grasp_gibbs_for_tuning(data = data,
                                    nsamples = pilot_nsamples,
                                    burnin = pilot_burnin,
                                    init = init,
                                    group_sizes = group_sizes,
                                    print_every = print_every,
                                    sigma_prop_a = s,
                                    sigma_prop_b = s,
                                    sigma_prop_a_global = s,
                                    sigma_prop_b_global = s,
                                    return_accept = TRUE)
    
    acc <- pilot$acceptance
    avg_a <- mean(acc$a_g, na.rm = TRUE)
    avg_b <- mean(acc$b_g, na.rm = TRUE)
    avg_a_global <- acc$a_global
    avg_b_global <- acc$b_global
    
    results[results$sigma == s, "acc_a"] <- avg_a
    results[results$sigma == s, "acc_b"] <- avg_b
    results[results$sigma == s, "acc_a_global"] <- avg_a_global
    results[results$sigma == s, "acc_b_global"] <- avg_b_global
    
    cat(sprintf("  Acceptance: a_g = %.3f, b_g = %.3f, a_global = %.3f, b_global = %.3f\n",
                avg_a, avg_b, avg_a_global, avg_b_global))
  }
  
  # Select best sigma for each type (closest to target_rate)
  best_idx_a <- which.min(abs(results$acc_a - target_rate))
  best_idx_b <- which.min(abs(results$acc_b - target_rate))
  best_idx_a_global <- which.min(abs(results$acc_a_global - target_rate))
  best_idx_b_global <- which.min(abs(results$acc_b_global - target_rate))
  
  best_sigma_a <- results$sigma[best_idx_a]
  best_sigma_b <- results$sigma[best_idx_b]
  best_sigma_a_global <- results$sigma[best_idx_a_global]
  best_sigma_b_global <- results$sigma[best_idx_b_global]
  
  cat("\n=============================================================\n")
  cat("Tuning complete.\n")
  cat(sprintf("Best sigma for a_g       : %.3f (acceptance %.3f)\n",
              best_sigma_a, results$acc_a[best_idx_a]))
  cat(sprintf("Best sigma for b_g       : %.3f (acceptance %.3f)\n",
              best_sigma_b, results$acc_b[best_idx_b]))
  cat(sprintf("Best sigma for a_global  : %.3f (acceptance %.3f)\n",
              best_sigma_a_global, results$acc_a_global[best_idx_a_global]))
  cat(sprintf("Best sigma for b_global  : %.3f (acceptance %.3f)\n",
              best_sigma_b_global, results$acc_b_global[best_idx_b_global]))
  cat("=============================================================\n")
  
  return(list(best_sigma_a = best_sigma_a,
              best_sigma_b = best_sigma_b,
              best_sigma_a_global = best_sigma_a_global,
              best_sigma_b_global = best_sigma_b_global,
              results = results))
}


# ==================================================
# Testing the function 
# ==================================================
# 
# Generate data (your own function)
dat <- generate_data(n = 500, setting = 1)
group_sizes <- dat$group_sizes
n_obs <- nrow(dat$data)
p <- ncol(dat$data) - 1
G <- length(group_sizes)

# Create initial vector (using your create_grasp_init)
init <- create_grasp_init(n_obs, p, G, group_sizes)


# Testing the grasp_gibbs_for_tuning

# grasp_gibbs_for_tuning(data = dat$data, nsamples = 100, burnin = 20, init= init,
#                        group_sizes = group_sizes,print_every = 50,
#                        sigma_prop_a = 0.2, sigma_prop_b = 0.2,
#                        sigma_prop_a_global = 0.2, sigma_prop_b_global = 0.2,
#                        return_accept = TRUE)



# Tune sigma values
tuned <- tune_sigma_grasp(data = dat$data,
                          init = init,
                          group_sizes = group_sizes,
                          sigma_candidates = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0,1.5,2,3,5),
                          target_rate = 0.3,
                          pilot_nsamples = 2000,
                          pilot_burnin = 500)
tuned


# ===================================================
# For every scenario
# ===================================================

settings = 1:4
tuned_sigmas = matrix(NA, nrow = 4, ncol = 4)
colnames(tuned_sigmas) = c("best_sigma_a", "best_sigma_b", "best_sigma_global_a",
                           "best_sigma_global_b")
rownames(tuned_sigmas)= paste0("Scenario :",1:4)

for(s in settings){
  cat(sprintf("==============================================="))
  cat(sprintf("\nScenario = %d ...\n", s))
  cat(sprintf("==============================================="))
  dat = generate_data(n = 500, setting = s)
  group_sizes = dat$group_sizes
  n_obs = nrow(dat$data)
  p = ncol(dat$data) - 1
  G = length(group_sizes)
  init = create_grasp_init(n_obs, p, G,group_sizes)
  
  tuned = tune_sigma_grasp(data = dat$data, init = init,group_sizes = group_sizes,
                           sigma_candidates = c(0.05,0.1,0.2,0.3,0.5,0.7,1.0,1.5,2,3,5),
                           target_rate = 0.3,
                           pilot_nsamples = 2000,
                           pilot_burnin = 500)
  tuned_sigmas[s,] = c(tuned$best_sigma_a,
                       tuned$best_sigma_b,
                       tuned$best_sigma_a_global,
                       tuned$best_sigma_b_global)
  
}

tuned_sigmas








# ===================================================
# For every scenario - For every correlation pair
# ===================================================


tune_scenario_grasp <- function(scenario, n = 500, n_replicates = 1,
                                sigma_candidates = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 3, 5),
                                target_rate = 0.3,
                                pilot_nsamples = 2000,
                                pilot_burnin = 500) {
  correlation_pairs <- list(
    c(rho_w = 0.8, rho_b = 0.2),
    c(rho_w = 0.5, rho_b = 0.1),
    c(rho_w = 0.3, rho_b = 0.01)
  )
  results_df <- data.frame(
    rho_w = sapply(correlation_pairs, function(x) x[1]),
    rho_b = sapply(correlation_pairs, function(x) x[2]),
    avg_sigma_a = NA,
    avg_sigma_b = NA,
    avg_sigma_a_global = NA,
    avg_sigma_b_global = NA
  )
  
  cat("\n===============================================\n")
  cat("Tuning GRASP Scenario", scenario, "\n")
  cat("===============================================\n")
  
  for (pair_idx in seq_along(correlation_pairs)) {
    rho_w <- correlation_pairs[[pair_idx]]["rho_w"]
    rho_b <- correlation_pairs[[pair_idx]]["rho_b"]
    
    cat(sprintf("\n--- Correlation: rho_w = %.2f, rho_b = %.2f ---\n", rho_w, rho_b))
    
    sum_sigma_a <- 0
    sum_sigma_b <- 0
    sum_sigma_a_global <- 0
    sum_sigma_b_global <- 0
    
    for (rep in 1:n_replicates) {
      cat(sprintf("  Replicate %d/%d\n", rep, n_replicates))
      
      # Generate data with current correlation
      dat <- generate_data(n = n, setting = scenario, rho_w = rho_w, rho_b = rho_b)
      group_sizes <- dat$group_sizes
      n_obs <- nrow(dat$data)
      p <- ncol(dat$data) - 1
      G <- length(group_sizes)
      init <- create_grasp_init(n_obs, p, G, group_sizes)   # Ensure this function exists
      
      tuned <- tune_sigma_grasp(
        data = dat$data,
        init = init,
        group_sizes = group_sizes,
        sigma_candidates = sigma_candidates,
        target_rate = target_rate,
        pilot_nsamples = pilot_nsamples,
        pilot_burnin = pilot_burnin,
        print_every = 1000
      )
      
      sum_sigma_a <- sum_sigma_a + tuned$best_sigma_a
      sum_sigma_b <- sum_sigma_b + tuned$best_sigma_b
      sum_sigma_a_global <- sum_sigma_a_global + tuned$best_sigma_a_global
      sum_sigma_b_global <- sum_sigma_b_global + tuned$best_sigma_b_global
    }
    
    avg_a <- sum_sigma_a / n_replicates
    avg_b <- sum_sigma_b / n_replicates
    avg_a_global <- sum_sigma_a_global / n_replicates
    avg_b_global <- sum_sigma_b_global / n_replicates
    
    results_df[pair_idx, "avg_sigma_a"] <- avg_a
    results_df[pair_idx, "avg_sigma_b"] <- avg_b
    results_df[pair_idx, "avg_sigma_a_global"] <- avg_a_global
    results_df[pair_idx, "avg_sigma_b_global"] <- avg_b_global
    
    cat(sprintf("  Average sigma: a = %.3f, b = %.3f, a_global = %.3f, b_global = %.3f\n",
                avg_a, avg_b, avg_a_global, avg_b_global))
  }
  
  cat("\n===============================================\n")
  cat("GRASP Scenario", scenario, "tuning complete.\n")
  print(results_df)
  cat("===============================================\n")
  
  return(results_df)
}



# Scenario 1
results_grasp_scenario1 <- tune_scenario_grasp(scenario = 1, n = 500, n_replicates = 1,
                                               sigma_candidates = c(0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 3, 5, 7, 10))



# Scenario 2
results_grasp_scenario2 <- tune_scenario_grasp(scenario = 2, n = 500, n_replicates = 1,
                                               sigma_candidates = c(0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 3, 5, 7, 10))

# Scenario 3 (n=300)
results_grasp_scenario3 <- tune_scenario_grasp(scenario = 3, n = 300, n_replicates = 1,
                                               sigma_candidates = c(0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 3, 5, 7, 10))

# Scenario 4 (n=200)
results_grasp_scenario4 <- tune_scenario_grasp(scenario = 4, n = 200, n_replicates = 1,
                                               sigma_candidates = c(0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 3, 5, 7, 10))



































