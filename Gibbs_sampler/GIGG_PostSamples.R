# ============================================================
# Helper to create initial vectors for GIGG (two types)
# ============================================================
create_gigg_init <- function(n, p, G, group_sizes, type = "zero") {
  # Order: w (n), beta (p), lambda_sq (p), delta_sq (G),
  #        tau_sq, v, a_sq (G), d (G), b_sq (G), e (G)
  if (type == "zero") {
    w <- rep(0.5, n)
    beta <- rep(0, p)
    lambda_sq <- rep(1, p)
    delta_sq <- rep(1, G)
    tau_sq <- 1
    v <- 1
    a_sq <- rep(1, G)
    d <- rep(1, G)
    b_sq <- rep(1, G)
    e <- rep(1, G)
  } else if (type == "overdispersed") {
    w <- rep(0.5, n)
    beta <- rnorm(p, 0, 1)
    lambda_sq <- runif(p, 0.5, 2)
    delta_sq <- runif(G, 0.5, 2)
    tau_sq <- runif(1, 0.5, 2)
    v <- runif(1, 0.5, 2)
    a_sq <- runif(G, 0.5, 2)
    d <- runif(G, 0.5, 2)
    b_sq <- runif(G, 0.5, 2)
    e <- runif(G, 0.5, 2)
  }
  return(c(w, beta, lambda_sq, delta_sq, tau_sq, v,
           a_sq, d, b_sq, e))
}

# ============================================================
# Function to run diagnostics for GIGG (similar to BGHS/GRASP)
# ============================================================
run_gigg_diagnostics <- function(scenario, rho_w, rho_b, n = 500,
                                 n_iter = 50000, burnin_keep = 20000,
                                 sigma_prop_a = 0.2, sigma_prop_b = 0.2,
                                 output_dir = "GIGG_diagnostics") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("\n==================================================\n")
  cat(sprintf("GIGG: Scenario %d, rho_w = %.1f, rho_b = %.2f\n", scenario, rho_w, rho_b))
  cat("==================================================\n")
  
  # Generate data
  dat <- generate_data(n = n, setting = scenario, rho_w = rho_w, rho_b = rho_b)
  group_sizes <- dat$group_sizes
  n_obs <- nrow(dat$data)
  p <- ncol(dat$data) - 1
  G <- length(group_sizes)
  
  # Initial values for two chains
  init1 <- create_gigg_init(n_obs, p, G, group_sizes, type = "zero")
  init2 <- create_gigg_init(n_obs, p, G, group_sizes, type = "overdispersed")
  
  # Run chains
  cat("  Running chain 1 (zero start)...\n")
  chain1 <- gigg_gibbs(data = dat$data, nsamples = n_iter, burnin = 0,
                       init = init1, group_sizes = group_sizes,
                       sigma_prop_a = sigma_prop_a,
                       sigma_prop_b = sigma_prop_b,
                       print_every = 5000)
  cat("  Running chain 2 (overdispersed start)...\n")
  chain2 <- gigg_gibbs(data = dat$data, nsamples = n_iter, burnin = 0,
                       init = init2, group_sizes = group_sizes,
                       sigma_prop_a = sigma_prop_a,
                       sigma_prop_b = sigma_prop_b,
                       print_every = 5000)
  
  # Beta columns
  beta_cols <- grep("^beta_", colnames(chain1))
  n_beta <- length(beta_cols)
  cat(sprintf("  Number of beta coefficients: %d\n", n_beta))
  
  # Create mcmc objects (only beta)
  mcmc1 <- as.mcmc(chain1[, beta_cols])
  mcmc2 <- as.mcmc(chain2[, beta_cols])
  mcmc_list <- mcmc.list(mcmc1, mcmc2)
  
  # Apply burn‑in
  mcmc_burn <- window(mcmc_list, start = burnin_keep + 1)
  
  # Univariate Gelman-Rubin
  rhat_univ <- gelman.diag(mcmc_burn, multivariate = FALSE)
  max_rhat <- max(rhat_univ$psrf[, "Point est."])
  
  # Multivariate Gelman-Rubin
  rhat_mult <- gelman.diag(mcmc_burn, multivariate = TRUE)
  mult_psrf <- rhat_mult$psrf
  cat(sprintf("  Multivariate PSRF: %.3f\n", mult_psrf))
  
  # Effective sample size (chain1 after burn‑in)
  ess <- effectiveSize(window(mcmc1, start = burnin_keep + 1))
  min_ess <- min(ess)
  cat(sprintf("  Min ESS: %.1f\n", min_ess))
  
  # ---------- Trace plots using coda::traceplot ----------
  plots_per_page <- 20
  n_pages <- ceiling(n_beta / plots_per_page)
  trace_file <- file.path(output_dir, sprintf("gigg_scen%d_rho%.1f_%.2f_trace_all_betas.pdf", scenario, rho_w, rho_b))
  pdf(trace_file, width = 12, height = 9)
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_beta)
    par(mfrow = c(4, 5), mar = c(2,2,2,1), oma = c(2,2,0,0))
    for (j in start_idx:end_idx) {
      single_list <- mcmc.list(mcmc1[, j - start_idx + 1],
                               mcmc2[, j - start_idx + 1])
      traceplot(single_list, main = colnames(chain1)[beta_cols[j]],
                xlab = "Iteration", ylab = "Value", col = c("red", "blue"))
    }
    mtext(sprintf("GIGG: Scenario %d, rho_w=%.1f, rho_b=%.2f - Trace plots (betas)", scenario, rho_w, rho_b),
          outer = TRUE, cex = 1.2, line = 0.5)
  }
  dev.off()
  cat(sprintf("  Trace plots saved to: %s\n", trace_file))
  
  # ---------- Autocorrelation plots ----------
  acf_file <- file.path(output_dir, sprintf("gigg_scen%d_rho%.1f_%.2f_acf_all_betas.pdf", scenario, rho_w, rho_b))
  pdf(acf_file, width = 12, height = 9)
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_beta)
    par(mfrow = c(4, 5), mar = c(3,3,2,1), oma = c(2,2,0,0))
    for (j in start_idx:end_idx) {
      autocorr.plot(window(mcmc1, start = burnin_keep + 1)[, j - start_idx + 1],
                    auto.layout = FALSE, main = colnames(chain1)[beta_cols[j]], lag.max = 50)
    }
    mtext(sprintf("GIGG: Scenario %d, rho_w=%.1f, rho_b=%.2f - ACF plots (betas, chain1)", scenario, rho_w, rho_b),
          outer = TRUE, cex = 1.2, line = 0.5)
  }
  dev.off()
  cat(sprintf("  ACF plots saved to: %s\n", acf_file))
  
  # Summary
  cat("\n--- Summary ---\n")
  cat(sprintf("GIGG Scenario %d, rho_w=%.1f, rho_b=%.2f\n", scenario, rho_w, rho_b))
  cat(sprintf("  Max univariate Rhat = %.3f %s\n", max_rhat, ifelse(max_rhat < 1.1, "(OK)", "(>1.1 - need more iterations)")))
  cat(sprintf("  Multivariate PSRF  = %.3f %s\n", mult_psrf, ifelse(mult_psrf < 1.1, "(OK)", "(>1.1 - need more iterations)")))
  cat(sprintf("  Min ESS            = %.1f %s\n", min_ess, ifelse(min_ess >= 400, "(OK)", "(<400 - need more samples)")))
  cat(sprintf("  Burn-in used = %d iterations\n", burnin_keep))
  cat(sprintf("  Total iterations per chain = %d\n", n_iter))
  cat("------------------------------------------------\n")
  
  return(list(scenario = scenario, rho_w = rho_w, rho_b = rho_b,
              max_rhat = max_rhat, mult_psrf = mult_psrf, min_ess = min_ess,
              n_beta = n_beta, burnin = burnin_keep, total_iter = n_iter))
}

# ============================================================
# Loop for GIGG diagnostics (using fixed tuned sigmas)
# ============================================================

# Define tuned sigmas for GIGG (replace with your actual values)
# You can make them scenario‑specific by using a vector/list per scenario.
# For simplicity, we use fixed values here.
sigma_prop_a_fixed <- 5   # your tuned sigma for a_g
sigma_prop_b_fixed <- 1.5   # your tuned sigma for b_g

# Define scenarios, sample sizes, and correlation pairs
scenario_n <- c(500, 500, 300, 200)
correlations <- list(
  c(rho_w = 0.8, rho_b = 0.2),
  c(rho_w = 0.5, rho_b = 0.1),
  c(rho_w = 0.3, rho_b = 0.01)
)

output_dir <- "GIGG_diagnostics"
all_results <- list()

# Start timing
total_start <- Sys.time()

for (scen in 1:4) {
  n_this <- scenario_n[scen]
  cat("===========================================================================\n")
  cat("\n\n========== Processing GIGG Scenario", scen, "with n =", n_this, "==========\n")
  cat("============================================================================\n")
  
  for (corr in correlations) {
    rho_w <- corr["rho_w"]
    rho_b <- corr["rho_b"]
    
    cat("\n--------------------------------------------------\n")
    cat(sprintf("Running for Scenario %d (n=%d), rho_w=%.1f, rho_b=%.2f\n", 
                scen, n_this, rho_w, rho_b))
    cat("--------------------------------------------------\n")
    
    result <- run_gigg_diagnostics(
      scenario = scen,
      rho_w = rho_w,
      rho_b = rho_b,
      n = n_this,
      n_iter = 20000,
      burnin_keep = 10000,
      sigma_prop_a = sigma_prop_a_fixed,
      sigma_prop_b = sigma_prop_b_fixed,
      output_dir = output_dir
    )
    
    all_results <- c(all_results, list(result))
  }
}

# End timing
total_end <- Sys.time()
elapsed_sec <- as.numeric(difftime(total_end, total_start, units = "secs"))
elapsed_min <- floor(elapsed_sec / 60)
elapsed_sec_rem <- round(elapsed_sec %% 60)
cat("\n========================================\n")
cat(sprintf("Total time for all GIGG diagnostics: %d minutes and %d seconds\n", 
            elapsed_min, elapsed_sec_rem))
cat("========================================\n")

# Convert to data frame and save summary
results_df <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    scenario = x$scenario,
    rho_w = x$rho_w,
    rho_b = x$rho_b,
    max_rhat = x$max_rhat,
    mult_psrf = x$mult_psrf,
    min_ess = x$min_ess,
    n_beta = x$n_beta,
    burnin_used = x$burnin,
    total_iter = x$total_iter
  )
}))

print(results_df)
write.csv(results_df, file.path(output_dir, "GIGG_diagnostics_summary.csv"), row.names = FALSE)

cat("\nAll GIGG diagnostics completed. Summary saved.\n")
