
# ================================================================
# GIGG Gibbs Sampler (with optional return of acceptance rates)
# ================================================================

gigg_gibbs_for_tuning <- function(data, nsamples = 1000, burnin = 1000,
                                  init, group_sizes, print_every = 100,
                                  sigma_prop_a = 0.2, sigma_prop_b = 0.2,
                                  return_accept = FALSE) {
  
  require(MASS); require(BayesLogit); require(MCMCpack); require(GIGrvg)
  
  n <- nrow(data)
  p <- ncol(data) - 1
  X <- as.matrix(data[, -1])
  y <- as.vector(data[, 1])
  G <- length(group_sizes)
  
  expected_len <- n + p + p + G + 1 + 1 + 4*G
  if(length(init) != expected_len) stop("init length mismatch")
  if(burnin < 0) stop("burnin non‑negative")
  if(burnin >= nsamples) stop("burnin < nsamples")
  if(sum(group_sizes) != p) stop("group_sizes sum != p")
  
  # Column names
  colname <- c(paste0("w_",1:n), paste0("beta_",1:p),
               unlist(lapply(1:G, function(g) paste0("lam_sq_", g, "_", 1:group_sizes[g]))),
               paste0("delta_sq_",1:G), "tau_sq", "v",
               paste0("a_sq_",1:G), paste0("d_",1:G),
               paste0("b_sq_",1:G), paste0("e_",1:G))
  
  sample <- matrix(NA, nrow = nsamples, ncol = length(colname))
  colnames(sample) <- colname
  
  # Extract initial values (order must match create_gigg_init)
  idx <- 0
  current_w <- init[(idx+1):(idx+n)]; idx <- idx + n
  current_beta <- init[(idx+1):(idx+p)]; idx <- idx + p
  current_lambda_sq <- init[(idx+1):(idx+p)]; idx <- idx + p
  current_delta_sq <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_tau_sq <- init[idx+1]; idx <- idx + 1
  current_v <- init[idx+1]; idx <- idx + 1
  current_a_sq <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_d <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_b_sq <- init[(idx+1):(idx+G)]; idx <- idx + G
  current_e <- init[(idx+1):(idx+G)]; idx <- idx + G
  
  idx_groups <- split(1:p, rep(1:G, times = group_sizes))
  k <- y - 0.5
  
  # Acceptance counters
  if (return_accept) {
    accept_a <- rep(0, G)
    accept_b <- rep(0, G)
    total_post_mh <- 0
  }
  
  # Initial row
  sample[1,] <- c(current_w, current_beta, current_lambda_sq,
                  current_delta_sq, current_tau_sq, current_v,
                  current_a_sq, current_d, current_b_sq, current_e)
  
  start_time <- Sys.time()
  cat("\n===============================================================\n")
  cat("           Starting GIGG Gibbs Sampler (Tuning)\n")
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
    current_w <- BayesLogit::rpg(n, 1, eta)
    
    # ---------- 2. Sample beta ----------
    D_vec <- current_tau_sq * current_delta_sq[rep(1:G, times = group_sizes)] * current_lambda_sq
    D_inv <- diag(1/D_vec)
    XtOmegaX <- crossprod(X, X * current_w)
    P <- XtOmegaX + D_inv
    R <- chol(P)
    mu <- backsolve(R, forwardsolve(t(R), crossprod(X, k)))
    current_beta <- mu + backsolve(R, rnorm(p))
    
    # ---------- 3. Sample lambda_sq ----------
    b_g_vec <- sqrt(current_b_sq[rep(1:G, times = group_sizes)])
    delta_g_vec <- rep(current_delta_sq, times = group_sizes)
    shape_lambda <- b_g_vec + 0.5
    scale_lambda <- 1 + (current_beta^2) / (2 * current_tau_sq * delta_g_vec)
    current_lambda_sq <- rinvgamma(p, shape = shape_lambda, scale = scale_lambda)
    
    # ---------- 4. Sample delta_sq (GIG) ----------
    new_delta_sq <- numeric(G)
    a_g_vec <- sqrt(current_a_sq)
    for(g in 1:G) {
      ind <- idx_groups[[g]]
      p_g <- group_sizes[g]
      sum_term <- sum(current_beta[ind]^2 / (current_tau_sq * current_lambda_sq[ind]))
      lambda_gig <- a_g_vec[g] - p_g/2
      chi_gig <- sum_term
      psi_gig <- 2.0
      new_delta_sq[g] <- rgig(1, lambda = lambda_gig, chi = chi_gig, psi = psi_gig)
    }
    current_delta_sq <- new_delta_sq
    
    # ---------- 5. Sample tau_sq and v ----------
    total_sum <- sum(current_beta^2 / (current_lambda_sq * rep(current_delta_sq, times = group_sizes)))
    shape_tau <- (p + 1) / 2
    scale_tau <- 0.5 * total_sum + 1 / current_v
    current_tau_sq <- rinvgamma(1, shape = shape_tau, scale = scale_tau)
    current_v <- rinvgamma(1, shape = 1, scale = 1 + 1 / current_tau_sq)
    
    # ---------- 6. Auxiliary variables ----------
    current_d <- rinvgamma(G, shape = 1, scale = 1 + 1 / current_a_sq)
    current_e <- rinvgamma(G, shape = 1, scale = 1 + 1 / current_b_sq)
    
    # ---------- 7. MH updates with acceptance counting ----------
    # a_sq
    for(g in 1:G) {
      old <- current_a_sq[g]
      new <- mh_update_a_g_sq_gigg(old, current_delta_sq[g], current_d[g], sigma_prop_a)
      if (return_accept && iter > burnin && new != old) accept_a[g] <- accept_a[g] + 1
      current_a_sq[g] <- new
    }
    # b_sq
    for(g in 1:G) {
      ind <- idx_groups[[g]]
      lambda_vals <- current_lambda_sq[ind]
      p_g <- group_sizes[g]
      old <- current_b_sq[g]
      new <- mh_update_b_g_sq_gigg(old, lambda_vals, current_e[g], p_g, sigma_prop_b)
      if (return_accept && iter > burnin && new != old) accept_b[g] <- accept_b[g] + 1
      current_b_sq[g] <- new
    }
    
    if (return_accept && iter > burnin) total_post_mh <- total_post_mh + 1
    
    # ---------- 8. Store samples ----------
    new_row <- c(current_w, current_beta, current_lambda_sq,
                 current_delta_sq, current_tau_sq, current_v,
                 current_a_sq, current_d, current_b_sq, current_e)
    sample[iter,] <- new_row
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
      a = accept_a / total_post_mh,
      b = accept_b / total_post_mh
    )
    return(list(posterior = posterior_sample, acceptance = acceptance))
  } else if (return_accept) {
    warning("No post‑burn‑in MH steps; acceptance rates NA")
    return(list(posterior = posterior_sample, acceptance = list(a = NA, b = NA)))
  } else {
    return(posterior_sample)
  }
}


# ============================================================
# Tune proposal sigmas for GRASP using a grid search
# ============================================================



tune_sigma_gigg = function(data, init, group_sizes,
                           sigma_candidates = c(0.05,0.1,0.2,0.3,0.5,
                                                0.7,1),
                           target_rate = 0.3,
                           pilot_nsamples = 500,
                           pilot_burnin = 300,
                           print_every = 200){
  n <- nrow(data)
  p <- ncol(data) - 1
  G <- length(group_sizes)
  
  results <- data.frame(sigma = sigma_candidates,
                        acc_a = NA, acc_b = NA)
  
  cat("\n=============================================================\n")
  cat("Tuning proposal sigma for GIGG\n")
  cat("=============================================================\n")
  cat(sprintf("Target acceptance rate: %.2f\n", target_rate))
  cat(sprintf("Pilot: %d iterations, burn‑in %d\n", pilot_nsamples, pilot_burnin))
  cat("Testing sigmas:", paste(sigma_candidates, collapse=", "), "\n")
  cat("=============================================================\n\n")
  
  for (s in sigma_candidates) {
    cat(sprintf("\nTesting sigma = %.3f ...\n", s))
    
    pilot <- gigg_gibbs_for_tuning(data = data,
                                   nsamples = pilot_nsamples,
                                   burnin = pilot_burnin,
                                   init = init,
                                   group_sizes = group_sizes,
                                   print_every = print_every,
                                   sigma_prop_a = s,
                                   sigma_prop_b = s,
                                   return_accept = TRUE)
    
    acc <- pilot$acceptance
    avg_a <- mean(acc$a, na.rm = TRUE)
    avg_b <- mean(acc$b, na.rm = TRUE)
    
    results[results$sigma == s, "acc_a"] <- avg_a
    results[results$sigma == s, "acc_b"] <- avg_b
    
    cat(sprintf("  Acceptance: a = %.3f, b = %.3f\n", avg_a, avg_b))
  }
  
  # Select best sigma for each type
  best_idx_a <- which.min(abs(results$acc_a - target_rate))
  best_idx_b <- which.min(abs(results$acc_b - target_rate))
  
  best_sigma_a <- results$sigma[best_idx_a]
  best_sigma_b <- results$sigma[best_idx_b]
  
  cat("\n=============================================================\n")
  cat("Tuning complete.\n")
  cat(sprintf("Best sigma for a_g       : %.3f (acceptance %.3f)\n",
              best_sigma_a, results$acc_a[best_idx_a]))
  cat(sprintf("Best sigma for b_g       : %.3f (acceptance %.3f)\n",
              best_sigma_b, results$acc_b[best_idx_b]))
  cat("=============================================================\n")
  
  return(list(best_sigma_a = best_sigma_a,
              best_sigma_b = best_sigma_b,
              results = results))
}



# ==================================================
# Testing the function 
# ==================================================
# 
# Generate data (your own function)
dat <- generate_data(n = 500, setting = 4)
group_sizes <- dat$group_sizes
n_obs <- nrow(dat$data)
p <- ncol(dat$data) - 1
G <- length(group_sizes)

# Create initial vector (using your create_grasp_init)
init <- create_gigg_init(n_obs, p, G, group_sizes)


# Testing the grasp_gibbs_for_tuning

# grasp_gibbs_for_tuning(data = dat$data, nsamples = 100, burnin = 20, init= init,
#                        group_sizes = group_sizes,print_every = 50,
#                        sigma_prop_a = 0.2, sigma_prop_b = 0.2,
#                        sigma_prop_a_global = 0.2, sigma_prop_b_global = 0.2,
#                        return_accept = TRUE)



# Tune sigma values
tuned <- tune_sigma_gigg(data = dat$data,
                          init = init,
                          group_sizes = group_sizes,
                          sigma_candidates = c(0.05, 0.1, 0.3, 0.5, 0.7, 1.0,1.5,2,3,5,6,7),
                          target_rate = 0.3,
                          pilot_nsamples = 2000,
                          pilot_burnin = 500)
tuned
























# ===================================================
# For every scenario
# ===================================================
# 
settings <- 1:4
tuned_sigmas_gigg <- matrix(NA, nrow = 4, ncol = 2)
colnames(tuned_sigmas_gigg) <- c("best_sigma_a", "best_sigma_b")
rownames(tuned_sigmas_gigg) <- paste0("Scenario_", settings)

# Optional: store elapsed times per scenario
elapsed_times <- numeric(length(settings))

for(s in settings) {
  cat("\n===============================================\n")
  cat(sprintf("Scenario = %d\n", s))
  cat("===============================================\n")

  start_time <- Sys.time()   # start timing

  dat <- generate_data(n = 500, setting = s)
  group_sizes <- dat$group_sizes
  n_obs <- nrow(dat$data)
  p <- ncol(dat$data) - 1
  G <- length(group_sizes)
  init <- create_gigg_init(n_obs, p, G, group_sizes)

  tuned <- tune_sigma_gigg(data = dat$data,
                           init = init,
                           group_sizes = group_sizes,
                           sigma_candidates = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2, 3, 5),
                           target_rate = 0.3,
                           pilot_nsamples = 2000,
                           pilot_burnin = 500)

  end_time <- Sys.time()     # end timing
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  elapsed_times[s] <- elapsed

  cat(sprintf("Time taken for scenario %d: %.1f seconds\n", s, elapsed))

  tuned_sigmas_gigg[s, ] <- c(tuned$best_sigma_a, tuned$best_sigma_b)
}

tuned_sigmas_gigg



















# ===================================================
# For every scenario - For every correlation pair
# ===================================================

tune_scenario_gigg = function(scenario, n = 500, n_replicates = 1,
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
    avg_sigma_b = NA
  )
  
  cat("\n===============================================\n")
  cat("Tuning Scenario", scenario, "\n")
  cat("===============================================\n")
  
  for (pair_idx in seq_along(correlation_pairs)) {
    rho_w <- correlation_pairs[[pair_idx]]["rho_w"]
    rho_b <- correlation_pairs[[pair_idx]]["rho_b"]
    
    cat(sprintf("\n--- Correlation: rho_w = %.2f, rho_b = %.2f ---\n", rho_w, rho_b))
    
    sum_sigma_a <- 0
    sum_sigma_b <- 0
    
    for (rep in 1:n_replicates) {
      cat(sprintf("  Replicate %d/%d\n", rep, n_replicates))
      
      # Generate data with current correlation
      dat <- generate_data(n = n, setting = scenario, rho_w = rho_w, rho_b = rho_b)
      group_sizes <- dat$group_sizes
      n_obs <- nrow(dat$data)
      p <- ncol(dat$data) - 1
      G <- length(group_sizes)
      init <- create_gigg_init(n_obs, p, G, group_sizes)
      
      tuned <- tune_sigma_gigg(
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
    }
    
    avg_a <- sum_sigma_a / n_replicates
    avg_b <- sum_sigma_b / n_replicates
    results_df[pair_idx, "avg_sigma_a"] <- avg_a
    results_df[pair_idx, "avg_sigma_b"] <- avg_b
    
    cat(sprintf("  Average sigma: a = %.3f, b = %.3f\n", avg_a, avg_b))
  }
  
  cat("\n===============================================\n")
  cat("Scenario", scenario, "tuning complete.\n")
  print(results_df)
  cat("===============================================\n")
  
  return(results_df)
  
}

results_gigg_scenario1 = tune_scenario_gigg(scenario = 1, n = 500, 
                                  n_replicates = 1,
                                  sigma_candidates = 
                                    c(0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 3, 5,7,10))
results_gigg_scenario1


results_gigg_scenario2 = tune_scenario_gigg(scenario = 2, n = 500, 
                                       n_replicates = 1,
                                       sigma_candidates = 
                                         c(0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 3, 5,7,10))
results_gigg_scenario2


results_gigg_scenario3 = tune_scenario_gigg(scenario = 3, n = 300, 
                                       n_replicates = 1,
                                       sigma_candidates = 
                                         c(0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 3, 5,7,10))
results_gigg_scenario3

results_gigg_scenario4 = tune_scenario_gigg(scenario = 4, n = 200, 
                                       n_replicates = 1,
                                       sigma_candidates = 
                                         c(0.05, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 3, 5,7,10))
results_gigg_scenario4
















