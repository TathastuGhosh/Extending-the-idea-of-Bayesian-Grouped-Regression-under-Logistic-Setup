

# ====================================================================
# -------------------- Gibbs sampling for BGHS -----------------------
# ====================================================================


bghs_gibbs = function(data, nsamples = 10000, burnin = 1000, 
                      init, group_sizes, print_every = 100){
  
  n = nrow(data)
  p = ncol(data) - 1
  X = as.matrix(data[,-1])
  y = as.vector(data[,1])
  G = length(group_sizes)
  expected_len = n + 3*p +2*G + 2
  
  # Few preliminary checks
  if(length(init) != expected_len) stop("init length does not match number of parameters.")
  if(burnin < 0) stop("burnin must be non-negative.")
  if(burnin >= nsamples) stop("burnin must be smaller nsamples.")
  if(sum(group_sizes) != p) stop("sum(group_sizes) must equal no. of covariates.")
  
  
  # Creating the sample dataframe
  colname = c(paste0("w_",1:n), paste0("beta_",1:p),
              unlist(Map(function(i,n) paste0("lam_sq_",i,1:n), 1:G, group_sizes)),
              unlist(Map(function(i,n) paste0("t_",i,1:n), 1:G, group_sizes)),
              paste0("delta_sq_",1:G),
              paste0("c_",1:G),
              "tau_sq","v"
  )
  sample = matrix(NA, nrow = nsamples, ncol = length(colname))
  colnames(sample) = colname
  
  # Tweaking initials
  current_w = init[1:n]
  current_beta = init[(n+1) :(n+p)]
  current_lambda_sq = init[(n+p+1) : (n+2*p)]
  current_t = init[(n+2*p+1) : (n+3*p)]
  current_delta_sq = init[(n+3*p+1) : (n+3*p+G)]
  current_c = init[(n+3*p+G+1) : (n+3*p+2*G)]
  current_tau_sq = init[n+3*p+2*G+1]
  current_v = init[n+3*p+2*G+2]
  
  # initial row
  sample[1,] = init 
  
  # Few predefined things
  k = y - 0.5
  idx = split(1:length(current_beta), rep(1:G, times = group_sizes)) # for generating delta square
  
  # Time related (Initial)
  start_time = Sys.time()
  cat("\n===============================================================\n")
  cat("                Starting BGHS Gibbs Sampler\n")
  cat("===============================================================\n")
  cat("\n\n")
  cat("----------------------------------------\n")
  cat("Summary of experimental setting\n")
  cat("----------------------------------------\n\n")
  cat(sprintf("Observations    : %d\n", n))
  cat(sprintf("Covariates      : %d\n", p))
  cat(sprintf("Number of groups: %d\n", G))
  cat(sprintf("Group sizes     : %s\n",paste(group_sizes,collapse = ", ")))
  
  cat("\n----------------------------------------\n\n")
  
  
  for(i in 2:nsamples){
    
    # progress printing
    if(i %% print_every == 0){
      elapsed = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      avg_time = elapsed/(i-1)
      eta_remaining = avg_time * (nsamples - i)
      pct = round(100*i/nsamples,1)
      current_time = format(Sys.time(),"%H:%M:%S")
      
      cat(sprintf(
        "Time: %s => Samples: [%d/%d] | Completed: %.1f%% | elapsed: %.1fs | ETA: %.1fs\n",
        current_time, i, nsamples, pct, elapsed, eta_remaining
      ))
    }
    
    
    # for generating w
    eta = as.vector(X %*% current_beta)
    current_w = BayesLogit::rpg(num= n, h = 1, z = eta)
    
    # for generating beta
    # XtOmegaX = crossprod(X, X*current_w)
    # D_inv = diag(1 / (rep(current_delta_sq,times = group_sizes) *current_lambda_sq *current_tau_sq))
    # Sigma_beta = solve(XtOmegaX + D_inv)
    # mu_beta = Sigma_beta %*% crossprod(X,k)
    # current_beta = MASS::mvrnorm(n = 1,
    #                              mu = as.vector(mu_beta),
    #                              Sigma = Sigma_beta)
    XtOmegaX = crossprod(X, X*current_w)
    D_inv = diag(1 / (rep(current_delta_sq,times = group_sizes) *current_lambda_sq *current_tau_sq))
    A = XtOmegaX + D_inv
    R = chol(A)
    Sigma_beta = chol2inv(R)
    mu_beta = Sigma_beta %*% crossprod(X,k)
    current_beta = MASS::mvrnorm(n = 1, mu = as.vector(mu_beta),
                                 Sigma = Sigma_beta)
    
    
    # for generating lambda square
    current_lambda_sq = MCMCpack::rinvgamma(p, shape = 1, scale = (current_beta^2)/(2*current_tau_sq*rep(current_delta_sq, group_sizes)) + (1/current_t))
    
    # for generating t 
    current_t = MCMCpack::rinvgamma(p, shape = 1, scale = 1 + (1/current_lambda_sq))
    
    # for generating delta square
    current_delta_sq = sapply(1:G,function(g){
      
      ind = idx[[g]]
      shape = (group_sizes[g] + 1)/2
      scale = (1/(2*current_tau_sq))*sum(current_beta[ind]^2 / current_lambda_sq[ind]) + 1/current_c[g]
      MCMCpack::rinvgamma(1, shape = shape, scale = scale)
    })
    
    # for generating c
    current_c = MCMCpack::rinvgamma(G, shape = 1, scale = 1 + (1/current_delta_sq))
    
    # for generating tau
    current_tau_sq = MCMCpack::rinvgamma(1, 
                                         shape = (p+1)/2, 
                                         scale = 
                                           0.5*sum(current_beta^2/
                                                     (current_lambda_sq * 
                                                        rep(current_delta_sq, 
                                                            times = group_sizes))) + 
                                           (1/current_v))
    
    
    # for generating v
    current_v = MCMCpack::rinvgamma(1,shape = 1, scale = 1 + 
                                      (1/current_tau_sq))
    
    new_row = c(current_w, current_beta, current_lambda_sq, current_t, current_delta_sq, current_c, current_tau_sq, current_v)
    sample[i,] = new_row # adding new row
  }
  
  sample = as.data.frame(sample)
  post_idx = seq(burnin + 1, nsamples, by = 1)
  posterior_sample = sample[post_idx,]
  rownames(posterior_sample) =seq_len(nrow(posterior_sample))
  
  # Time related (End)
  elapsed =as.numeric(difftime(Sys.time(),start_time,units = "secs"))
  cat("\n----------------------------------------\n")
  cat("Run Summary\n")
  cat("----------------------------------------\n")
  cat(sprintf("Total samples : %d\n",nsamples))
  cat(sprintf("Burn-in       : %d\n",burnin))
  cat(sprintf("Posterior kept: %d\n",nrow(posterior_sample)))
  elapsed_min = floor(elapsed / 60)
  elapsed_sec = round(elapsed %% 60)
  cat(sprintf("Elapsed time  : %.1f seconds (%dm %02ds)\n",elapsed, elapsed_min, elapsed_sec))
  cat("----------------------------------------\n\n")
  return(posterior_sample)
}





# ------------------------------
# checking
# -------------------------------


# result = generate_data(n = 500, rho_w = 0.8, rho_b = 0.2, setting = 2)
# data = result$data
# group_sizes = result$group_sizes
# tot_param = nrow(data) + 3*(ncol(data)-1) + 2*length(group_sizes) + 2
# init = c(rep(0.25,nrow(data)), rep(0,sum(group_sizes)),rep(1, tot_param-(nrow(data)+sum(group_sizes))))
# 
# post_sample_BGHS = bghs_gibbs(data = data,nsamples = 80000, burnin = 15000,
#                          init = init, group_sizes = group_sizes,
#                          print_every = 1000)
# 
# 
# 
# sample2 = bghs_gibbs(data = data,nsamples = 50000, burnin = 15000,
#                      init = init, group_sizes = group_sizes,
#                      print_every = 5000)
# 
# 
# sample3 = bghs_gibbs(data = data,nsamples = 500, burnin = 50,
#                      init = init, group_sizes = group_sizes,
#                      print_every = 50)


# -----------------------


res = generate_data(n = 500, rho_w = 0.8, rho_b = 0.2, setting = 1)
data = res$data
group_sizes = res$group_sizes
n = nrow(data)
p = ncol(data) - 1
tot_param = nrow(data) + 3*(ncol(data)-1) + 2*length(group_sizes) + 2

# Creating initials

init1 = c(rep(0.25,nrow(data)), rep(0,sum(group_sizes)),rep(1, tot_param-(nrow(data)+sum(group_sizes))))

init2 = init1
init2[1:n] = runif(n,0.1,1)
init2[(n+1) : (n+p)] = rnorm(p, mean = 0, sd = 1)
init2[(n+p+1):length(init1)] = runif(length(init1)-n-p,0.5,3)

init3 = init1
init3[1:n] = runif(n,1,5)
init3[(n+1):(n+p)] = rnorm(p, mean = 0, sd = 1)
init3[(n+p+1):length(init1)] = runif(length(init1)-n-p, 3,10)

# Drawing samples

sample1 = bghs_gibbs(data = data, nsamples = 100000, burnin = 0, init = init1,
                     print_every = 10000, group_sizes = group_sizes)
sample2 = bghs_gibbs(data = data, nsamples = 100000, burnin = 0, init = init2,
                     print_every = 10000, group_sizes = group_sizes)
sample3 = bghs_gibbs(data = data, nsamples = 100000, burnin = 0, init = init3,
                     print_every = 10000, group_sizes = group_sizes)

# coda thing

library(coda)

choose_burnin_GR = function(sample1, sample2, sample3,
                            threshold = 1.1,
                            by= 500,
                            max_burnin = 50000){
  chain1 = mcmc(sample1)
  chain2 = mcmc(sample2)
  chain3 = mcmc(sample3)
  
  # beta_i, lambda_sq_i, delta_sq_i, tau_sq, v
  param_names = colnames(sample1)
  monitored_params = c(
    grep("^beta_", param_names, value = TRUE),
    #grep("^lam_sq_", param_names, value = TRUE)
    grep("^delta_sq_", param_names, value = TRUE),
    "tau_sq"
  )
  burnin_grid = seq(0, max_burnin, by = by)
  results = data.frame(burnin = burnin_grid,
                       max_rhat = NA,
                       mean_rhat = NA,
                       all_converged = NA)
  chosen_burnin = NA
  final_rhat = NULL
  
  for(i in seq_along(burnin_grid)){
    b = burnin_grid[i]
    mcmc_post = mcmc.list(
      window(chain1[,monitored_params], start = b+1),
      window(chain2[,monitored_params], start = b+1),
      window(chain3[,monitored_params], start = b+1)
    )
    rhat = gelman.diag(mcmc_post, multivariate = FALSE,
                       autoburnin = FALSE)$psrf[,1]
    rhat = rhat[is.finite(rhat)]
    results$max_rhat[i] = max(rhat)
    results$mean_rhat[i] = mean(rhat)
    
    # TRUE only if EVERY parameter satisifies GR criterion
    results$all_converged[i] = all(rhat < threshold)
    cat("Burn-in =",b,"| Max Rhat =",round(max(rhat),4),"\n")
    
    if(is.na(chosen_burnin) && all(rhat < threshold)){
      chosen_burnin = b
      final_rhat = rhat
    }
  }
  return(list(diagnostics = results,
              monitored_parameters = monitored_params,
              rhat = final_rhat))
}

burnin_out = choose_burnin_GR(
  sample1, sample2, sample3,
  threshold = 1.1, # 1.08
  by = 5000,
  max_burnin = 90000
)
burnin_out$monitored_parameters
burnin_out$diagnostics





burnin = 15000 # setting burnin


# Samples after burnin
mcmc_post = mcmc.list(
  window(mcmc(sample1), start = burnin + 1),
  window(mcmc(sample2), start = burnin + 1),
  window(mcmc(sample3), start = burnin + 1)
)


# Extracting the important index
param_names = colnames(sample1)

monitor_idx = c(
  grep("^beta_", param_names),
  grep("^lam_sq_", param_names),
  grep("^delta_sq_", param_names),
  which(param_names == "tau_sq"),
  which(param_names == "v")
)

# Gelman-Rubin for the above index
gr_main =
  gelman.diag(
    mcmc_post[, monitor_idx],
    autoburnin = FALSE,
    multivariate = TRUE
  )$psrf[,1]

summary(gr_main)
sort(gr_main, decreasing = TRUE)[1:20]


# ESS for the above index
ess_main =
  effectiveSize(
    mcmc_post[, monitor_idx]
  )

summary(ess_main)
min(ess_main)
median(ess_main)


# Observing Gelman-Rubin and ESS for the lambda
lambda_idx =
  grep("^lam_sq_", param_names)

lambda_rhat =
  gelman.diag(
    mcmc_post[, lambda_idx],
    autoburnin = FALSE,
    multivariate = FALSE
  )$psrf[,1]
summary(gr_main)

ess_lambda =
  effectiveSize(
    mcmc_post[, lambda_idx]
  )

summary(ess_lambda)

min(ess_lambda)

median(ess_lambda)

worst_lambda =
  names(
    sort(lambda_rhat,
         decreasing = TRUE)[1:20]
  )

worst_lambda


par(mfrow = c(3,3))

traceplot(
  mcmc_post[, worst_lambda]
)


# acf plots for betas
param_names = colnames(sample1)

beta_names =
  grep("^beta_",
       param_names,
       value = TRUE)

par(mfrow = c(3,3))

for(b in beta_names[1:9]){
  
  acf(
    as.numeric(mcmc_post[[1]][, b]),
    lag.max = 100,
    main = b
  )
}

par(mfrow = c(1,1))

beta_idx =
  grep("^beta_", colnames(sample1))

ess_beta =
  effectiveSize(
    mcmc_post[, beta_idx]
  )
summary(ess_beta)

min(ess_beta)

median(ess_beta)


# traceplots

set.seed(123)

# beta parameter names
beta_names =
  grep("^beta_",
       colnames(sample1),
       value = TRUE)

# randomly choose 9 betas
beta_subset =
  sample(beta_names, 9)

par(mfrow = c(3,3))

traceplot(
  mcmc_post[, beta_subset]
)

par(mfrow = c(1,1))



















