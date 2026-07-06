# =====================================================================
# =====================================================================
# 
# ========================= GRASP SCENARIO ============================
#  
# =====================================================================
# =====================================================================

# USER INPUT (Change these for each run)
# -------------------------------------------------
scenario <- 1                         # 1,2,3,4 
rho_w <- 0.5                          # Within group corr : 0.8, 0.5, 0.3
rho_b <- 0.3                          # between group corr : 0.2, 0.1, 0.01
n <- c(500,500,300,200)[scenario]     # sample size
n_iter <- 20000                       # total mcmc per chain
burnin <- 10000                       # burnin samples
n_reps <- 10                          # replications

# Tuned proposal sigmas for GRASP (replace with your actual tuned values)
sigma_prop_a <- 1.5
sigma_prop_b <- 1.5
sigma_prop_a_global <- 3
sigma_prop_b_global <- 1.5

# output file (different from BGHS to avoid mixing)
results_file <- "GRASP_results.csv"

# No further changes below this line
# ------------------------------------

library(MASS)
library(BayesLogit)
library(MCMCpack)

# Helper: create initial vector for GRASP (zero start)
create_grasp_init <- function(n, p, G, group_sizes, type = "zero") {
  # Order: w (n), beta (p), lambda_sq (p), t (p), delta_sq (G), c (G),
  #        tau_sq, v, a_sq (G), d (G), b_sq (G), e (G),
  #        a_sq_global, d_global, b_sq_global, e_global
  if (type == "zero") {
    w <- rep(0.5, n)
    beta <- rep(0, p)
    lambda_sq <- rep(1, p)
    t <- rep(1, p)
    delta_sq <- rep(1, G)
    c <- rep(1, G)
    tau_sq <- 1
    v <- 1
    a_sq <- rep(1, G)
    d <- rep(1, G)
    b_sq <- rep(1, G)
    e <- rep(1, G)
    a_global <- 1
    d_global <- 1
    b_global <- 1
    e_global <- 1
  } else {
    # Overdispersed (optional – we keep zero for simplicity)
    w <- rep(0.5, n)
    beta <- rnorm(p, 0, 1)
    lambda_sq <- runif(p, 0.5, 2)
    t <- runif(p, 0.5, 2)
    delta_sq <- runif(G, 0.5, 2)
    c <- runif(G, 0.5, 2)
    tau_sq <- runif(1, 0.5, 2)
    v <- runif(1, 0.5, 2)
    a_sq <- runif(G, 0.5, 2)
    d <- runif(G, 0.5, 2)
    b_sq <- runif(G, 0.5, 2)
    e <- runif(G, 0.5, 2)
    a_global <- runif(1, 0.5, 2)
    d_global <- runif(1, 0.5, 2)
    b_global <- runif(1, 0.5, 2)
    e_global <- runif(1, 0.5, 2)
  }
  return(c(w, beta, lambda_sq, t, delta_sq, c, tau_sq, v,
           a_sq, d, b_sq, e, a_global, d_global, b_global, e_global))
}

# --------------------------------------------------
# Start Timing
# --------------------------------------------------
total_start <- Sys.time()

# --------------------------------------------------
# Simulation Loop
# --------------------------------------------------
cat("\n============================================================\n")
cat(sprintf("GRASP: Scenario %d, rho_w=%.1f, rho_b=%.2f, n=%d\n", scenario, rho_w, rho_b, n))
cat(sprintf("MCMC: %d total, %d burn‑in, %d replications\n", n_iter, burnin, n_reps))
cat("============================================================\n")

# Storage vectors for MSE
mse_null_vec <- numeric(n_reps)
mse_nonnull_vec <- numeric(n_reps)
mse_total_vec <- numeric(n_reps)

for(rep in 1:n_reps){
  
  cat(sprintf("-------------------------------------------------------------------------------------\n"))
  cat(sprintf("-------------------------------------------------------------------------------------\n"))
  cat(sprintf(" REPLICATION : %d/%d\n", rep, n_reps))
  cat(sprintf("-------------------------------------------------------------------------------------\n"))
  cat(sprintf("-------------------------------------------------------------------------------------\n"))
  
  # 1. Generate a new dataset (uses your existing generate_data)
  dat <- generate_data(n = n, rho_w = rho_w, rho_b = rho_b, setting = scenario)
  group_sizes <- dat$group_sizes
  true_beta <- dat$beta
  n_obs <- nrow(dat$data)
  p <- ncol(dat$data) - 1
  G <- length(group_sizes)
  
  # 2. Initial values (zero start)
  init <- create_grasp_init(n_obs, p, G, group_sizes, type = "zero")
  
  # 3. Run GRASP Gibbs sampler
  posterior <- grasp_gibbs(data = dat$data,
                           nsamples = n_iter,
                           burnin = burnin,
                           init = init,
                           group_sizes = group_sizes,
                           sigma_prop_a = sigma_prop_a,
                           sigma_prop_b = sigma_prop_b,
                           sigma_prop_a_global = sigma_prop_a_global,
                           sigma_prop_b_global = sigma_prop_b_global,
                           print_every = 1000)
  
  # 4. Extract posterior mean of regression coefficients
  beta_cols <- grep("^beta_", colnames(posterior))
  beta_hat <- colMeans(posterior[, beta_cols])
  
  # 5. Identify null and non‑null indices from the true beta
  idx_null <- which(true_beta == 0)
  idx_nonnull <- which(true_beta != 0)
  
  # 6. Compute MSE for this replication
  mse_null_vec[rep] <- mean((beta_hat[idx_null] - 0)^2)
  mse_nonnull_vec[rep] <- mean((beta_hat[idx_nonnull] - true_beta[idx_nonnull])^2)
  mse_total_vec[rep] <- mean((beta_hat - true_beta)^2)
}

# ----------------------------------------------------------------------
# End Timing
# ----------------------------------------------------------------------
total_end <- Sys.time()
elapsed_sec <- as.numeric(difftime(total_end, total_start, units = "secs"))
elapsed_min <- floor(elapsed_sec / 60)
elapsed_sec_rem <- round(elapsed_sec %% 60)

# ----------------------------------------------------------------------
# Aggregate Results
# ----------------------------------------------------------------------
mean_null <- mean(mse_null_vec)
se_null <- sd(mse_null_vec)
mean_nonnull <- mean(mse_nonnull_vec)
se_nonnull <- sd(mse_nonnull_vec)
mean_total <- mean(mse_total_vec)
se_total <- sd(mse_total_vec)

# Print to console
cat("\n------------------------------------------------------------\n")
cat(sprintf("Results for scenario %d, rho_w = %.1f, rho_b = %.2f\n", scenario, rho_w, rho_b))
cat(sprintf("MSE null     = %.6f (%.6f)\n", mean_null, se_null))
cat(sprintf("MSE non‑null = %.6f (%.6f)\n", mean_nonnull, se_nonnull))
cat(sprintf("MSE total    = %.6f (%.6f)\n", mean_total, se_total))
cat(sprintf("\nTotal elapsed time: %d minutes and %d seconds\n", elapsed_min, elapsed_sec_rem))
cat("------------------------------------------------------------\n")

# -------------------------------------------------------------------
# Append results to CSV file
# -------------------------------------------------------------------
new_row <- data.frame(
  method = "GRASP",
  scenario = scenario,
  rho_w = rho_w,
  rho_b = rho_b,
  n = n,
  n_iter = n_iter,
  burnin = burnin,
  n_reps = n_reps,
  MSE_null = sprintf("%.6f (%.6f)", mean_null, se_null),
  MSE_nonnull = sprintf("%.6f (%.6f)", mean_nonnull, se_nonnull),
  MSE_total = sprintf("%.6f (%.6f)", mean_total, se_total)
)

if (!file.exists(results_file)) {
  write.csv(new_row, results_file, row.names = FALSE)
  cat(sprintf("\nCreated new file: %s\n", results_file))
} else {
  existing <- read.csv(results_file)
  # Avoid duplicate entries for the same (method, scenario, rho_w, rho_b)
  dup <- duplicated(rbind(existing[,1:4], new_row[,1:4]))[nrow(existing)+1]
  if (!dup) {
    combined <- rbind(existing, new_row)
    write.csv(combined, results_file, row.names = FALSE)
    cat(sprintf("\nAppended results to: %s\n", results_file))
  } else {
    cat(sprintf("\nCombination already exists in %s – not appended.\n", results_file))
  }
}

