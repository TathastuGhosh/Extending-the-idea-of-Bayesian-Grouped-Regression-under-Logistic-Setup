
# ============================================================
# Load required libraries
# ============================================================
library(MASS)
library(BayesLogit)
library(MCMCpack)
library(coda)

# ============================================================
# Helper: Create three different initial vectors for BGHS
# ============================================================
create_bghs_init <- function(n, p, G, group_sizes, type = "zero") {
  # Order: w (n), beta (p), lambda_sq (p), t (p), delta_sq (G), c (G), tau_sq, v
  if (type == "zero") {
    w <- rep(0.5, n)
    beta <- rep(0, p)
    lambda_sq <- rep(1, p)
    t <- rep(1, p)
    delta_sq <- rep(1, G)
    c <- rep(1, G)
    tau_sq <- 1
    v <- 1
  } else if (type == "random_small") {
    w <- rep(0.5, n)
    beta <- rnorm(p, 0, 0.1)
    lambda_sq <- rep(1, p)
    t <- rep(1, p)
    delta_sq <- rep(1, G)
    c <- rep(1, G)
    tau_sq <- 1
    v <- 1
  } else if (type == "overdispersed") {
    w <- rep(0.5, n)
    beta <- rnorm(p, 0, 1)
    lambda_sq <- runif(p, 0.5, 2)
    t <- runif(p, 0.5, 2)
    delta_sq <- runif(G, 0.5, 2)
    c <- runif(G, 0.5, 2)
    tau_sq <- runif(1, 0.5, 2)
    v <- runif(1, 0.5, 2)
  }
  return(c(w, beta, lambda_sq, t, delta_sq, c, tau_sq, v))
}

# ============================================================
# Function to run diagnostics for a single combination
# ============================================================

run_bghs_diagnostics <- function(scenario, rho_w, rho_b, n = 500,
                                 n_iter = 80000, burnin_keep = 20000,
                                 output_dir = "BGHS_diagnostics") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("\n==================================================\n")
  cat(sprintf("Scenario %d, rho_w = %.1f, rho_b = %.2f\n", scenario, rho_w, rho_b))
  cat("==================================================\n")
  
  # Generate data
  dat <- generate_data(n = n, setting = scenario, rho_w = rho_w, rho_b = rho_b)
  group_sizes <- dat$group_sizes
  n_obs <- nrow(dat$data)
  p <- ncol(dat$data) - 1
  G <- length(group_sizes)
  
  # Initial values
  init1 <- create_bghs_init(n_obs, p, G, group_sizes, type = "zero")
  init2 <- create_bghs_init(n_obs, p, G, group_sizes, type = "overdispersed")
  
  # Run chains
  cat("  Running chain 1 (zero start)...\n")
  chain1 <- bghs_gibbs(data = dat$data, nsamples = n_iter, burnin = 0,
                       init = init1, group_sizes = group_sizes, print_every = 2000)
  cat("  Running chain 2 (overdispersed start)...\n")
  chain2 <- bghs_gibbs(data = dat$data, nsamples = n_iter, burnin = 0,
                       init = init2, group_sizes = group_sizes, print_every = 2000)
  
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
  trace_file <- file.path(output_dir, sprintf("scen%d_rho%.1f_%.2f_trace_all_betas.pdf", scenario, rho_w, rho_b))
  pdf(trace_file, width = 12, height = 9)
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_beta)
    # Set layout: 4 rows, 5 columns
    par(mfrow = c(4, 5), mar = c(2,2,2,1), oma = c(2,2,0,0))
    for (j in start_idx:end_idx) {
      # Create an mcmc.list for this single coefficient (all three chains)
      single_list <- mcmc.list(mcmc1[, j - start_idx + 1],
                               mcmc2[, j - start_idx + 1])
      # Use coda::traceplot (overlays chains with different colours)
      traceplot(single_list, main = colnames(chain1)[beta_cols[j]],
                xlab = "Iteration", ylab = "Value", col = c("red", "blue"))
    }
    mtext(sprintf("Scenario %d, rho_w=%.1f, rho_b=%.2f - Trace plots (betas)", scenario, rho_w, rho_b),
          outer = TRUE, cex = 1.2, line = 0.5)
  }
  dev.off()
  cat(sprintf("  Trace plots saved to: %s\n", trace_file))
  
  # ---------- Autocorrelation plots using coda::autocorr.plot ----------
  acf_file <- file.path(output_dir, sprintf("scen%d_rho%.1f_%.2f_acf_all_betas.pdf", scenario, rho_w, rho_b))
  pdf(acf_file, width = 12, height = 9)
  for (page in 1:n_pages) {
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_beta)
    par(mfrow = c(4, 5), mar = c(3,3,2,1), oma = c(2,2,0,0))
    for (j in start_idx:end_idx) {
      # Use chain1 after burn‑in
      autocorr.plot(window(mcmc1, start = burnin_keep + 1)[, j - start_idx + 1],
                    auto.layout = FALSE, main = colnames(chain1)[beta_cols[j]], lag.max = 50)
    }
    mtext(sprintf("Scenario %d, rho_w=%.1f, rho_b=%.2f - ACF plots (betas, chain1)", scenario, rho_w, rho_b),
          outer = TRUE, cex = 1.2, line = 0.5)
  }
  dev.off()
  cat(sprintf("  ACF plots saved to: %s\n", acf_file))
  
  # Summary
  cat("\n--- Summary ---\n")
  cat(sprintf("Scenario %d, rho_w=%.1f, rho_b=%.2f\n", scenario, rho_w, rho_b))
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







# Testing
# ---------------------------------


run_bghs_diagnostics(scenario = 1, rho_w = 0.5, rho_b = 0.1,
                     n = 500, n_iter = 100, burnin_keep = 20,
                     output_dir = "BGHS_test")






result = run_bghs_diagnostics(scenario = 1,
                              rho_w = 0.8,
                              rho_b = 0.2,
                              n = 500,
                              n_iter = 50000,
                              burnin_keep = 20000,
                              output_dir = "BGHS_diagnostic_1")





# Ultimatum
# ---------------------------------

# Define scenarios, their sample sizes, and correlation pairs
scenario_n <- c(500, 500, 300, 200)   # n for scenarios 1,2,3,4
correlations <- list(
  c(rho_w = 0.8, rho_b = 0.2),
  c(rho_w = 0.5, rho_b = 0.1),
  c(rho_w = 0.3, rho_b = 0.01)
)

output_dir <- "BGHS_diagnostics"
all_results <- list()

# start_timing
total_start = Sys.time()


# Loop over scenarios 1..4
for (scen in 1:4) {
  n_this <- scenario_n[scen]
  cat("===========================================================================\n")
  cat("\n\n========== Processing Scenario", scen, "with n =", n_this, "==========\n")
  cat("============================================================================\n")
  
  for (corr in correlations) {
    rho_w <- corr["rho_w"]
    rho_b <- corr["rho_b"]
    
    cat("\n--------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat(sprintf("Running for Scenario %d (n=%d), rho_w=%.1f, rho_b=%.2f\n", 
                scen, n_this, rho_w, rho_b))
    cat("-----------------------------------------------------\n")
    cat("--------------------------------------------------\n")
    
    result <- run_bghs_diagnostics(
      scenario = scen,
      rho_w = rho_w,
      rho_b = rho_b,
      n = n_this,               # scenario‑specific sample size
      n_iter = 20000,
      burnin_keep = 10000,
      output_dir = output_dir
    )
    
    all_results <- c(all_results, list(result))
  }
}
# End timing
total_end = Sys.time()
elapsed_sec <- as.numeric(difftime(total_end, total_start, units = "secs"))
elapsed_min <- floor(elapsed_sec / 60)
elapsed_sec_rem <- round(elapsed_sec %% 60)
cat("\n========================================\n")
cat(sprintf("Total time for all scenarios: %d minutes and %d seconds\n", 
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
write.csv(results_df, file.path(output_dir, "BGHS_diagnostics_summary.csv"), row.names = FALSE)

cat("\nAll diagnostics completed. Summary saved.\n")




















