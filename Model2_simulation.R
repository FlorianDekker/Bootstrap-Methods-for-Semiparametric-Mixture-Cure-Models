# Semiparametric Mixture Cure Model Simulation Study
# 
# This script simulates survival data under a semiparametric mixture cure model with a Weibull PH latency function,
# applies four different bootstrap methods (naive, stratified, conditional1, conditional2),
# and evaluates confidence interval coverage and interval length.

## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------------------------------------------------------------------
generate_weibull_PH <- function(x, mu = 1.5, rho = 1.75, beta = c(1, -0.1, 0.8), tau0 = 4) {
  # n: number of observations (uncured)
  # x: vector of covariates for latency (length = n)
  
  n <- nrow(x)
  u <- runif(n)
  T <- (-log(u) / (mu * exp(x %*% beta)))^(1 / rho)
  
  # Apply truncation at tau0
  T_trunc <- pmin(T, tau0)
  
  return(T_trunc)
}

## ------------------------------------------------------------------------------------------------------------------------------------
simulate_data <- function(n, beta = c(1, -0.1, 0.8), gamma = c(-0.8, 1.3, 1.5, -0.2), lambda_C = 1.3, tau = 7) {
  Z_1 <- matrix(rnorm(n, 0, 2), ncol = 1)
  Z_2 <- matrix(rbinom(n, size = 1, prob = 0.6), ncol = 1)
  Z_3 <- matrix(rbinom(n, size = 1, prob = 0.4), ncol = 1)
  Z <- cbind(rep(1,n), Z_1, Z_2, Z_3)
  X_1 <- Z_1
  X_2 <- matrix(runif(n, -3, 3), ncol = 1)
  X_3 <- Z_2
  X <- cbind(X_1, X_2, X_3)
  p_cure <- exp(Z %*% gamma) / (1 + exp(Z %*% gamma))
  cure_status <- rbinom(n, 1, p_cure)  # 1 = uncured
  T <- rep(Inf, n)
  uncured_idx <- which(cure_status == 1)
  T[uncured_idx] <- generate_weibull_PH(X[uncured_idx,])
  
  censor_times <- rexp(n, rate = lambda_C)
  censor_times_trunc <- pmin(censor_times, tau)
  censored <- as.numeric(censor_times_trunc >= T)
  
  time <- pmin(T, censor_times_trunc)
  
  data.frame(time = time, censored = censored, cure_status = cure_status, X = X, Z = Z)
}


## ------------------------------------------------------------------------------------------------------------------------------------
# Master bootstrap wrapper
do_bootstrap <- function(fit, method = "naive", nboot = 500, eps = 1e-7, link = "logit", emmax = 50) {
  if (method == "naive") {
    return(naive_bootstrap(fit, nboot, eps, link, emmax))
  }
  if (method == "naive stratified") {
    return(naive_bootstrap_stratified(fit, nboot, eps, link, emmax))
  }
  if (method == "conditional1") {
    return(conditional_bootstrap_1(fit, nboot, eps, link, emmax))
  }
  if (method == "conditional2") {
    return(conditional_bootstrap_2(fit, nboot, eps, link, emmax))
  }
  stop("Unknown bootstrap method: ", method)
}


## ------------------------------------------------------------------------------------------------------------------------------------
naive_bootstrap <- function(fit, nboot = 500, eps = 1e-7, link = "logit", emmax = 50) {
  n <- nrow(fit$X)
  nb <- length(fit$bnm)
  nbeta <- length(fit$betanm)
  censoring_rates <- numeric(nboot)
  percentage_plateau <- numeric(nboot)
  max_event_times <- numeric(nboot)
  
  tempdata <- cbind(fit$Time, fit$Status, fit$X, fit$Z)
  colnames(tempdata)[1:2] <- c("Time", "Status")  
  
  b_boot <- matrix(0, nboot, nb)
  beta_boot <- matrix(0, nboot, nbeta)
  
  i <- 1
  while (i <= nboot) {
    # Non-stratified resampling: sample n rows from the full dataset
    bootdata <- tempdata[sample(1:n, n, replace = TRUE), ]
    
    bootZ <- bootdata[, fit$bnm]
    bootX <- as.matrix(cbind(rep(1, n), bootdata[, fit$betanm]))
    
    bootfit <- em(bootdata[, "Time"], bootdata[, "Status"], bootX, bootZ, NULL, fit$b, fit$beta, link, emmax, eps)
    
    if (bootfit$tau < eps) {
      b_boot[i, ] <- bootfit$b
      beta_boot[i, ] <- bootfit$latencyfit
      Delta_star <- bootdata[, "Status"]
      Y_star <- bootdata[, "Time"]
      censoring_rates[i] <- 1 - mean(Delta_star)
      
      last_event_time <- max(Y_star[Delta_star == 1])
      max_event_times[i] <- last_event_time
      censored_after <- sum(Y_star[Delta_star == 0] > last_event_time)
      percentage_plateau[i] <- censored_after / n
      
      i <- i + 1
    }
  }
  
  b_var <- apply(b_boot, 2, var)
  beta_var <- apply(beta_boot, 2, var)
  b_sd <- sqrt(b_var)
  beta_sd <- sqrt(beta_var)
  
  list(
    b_boot = b_boot,
    beta_boot = beta_boot,
    b_sd = b_sd,
    beta_sd = beta_sd,
    censoring_rates = censoring_rates,
    percentage_plateau = percentage_plateau,
    max_event_times = max_event_times,
    error = bootfit$error
  )
}

## ------------------------------------------------------------------------------------------------------------------------------------
# Naive stratified bootstrap function
naive_bootstrap_stratified <- function(fit, nboot = 500, eps = 1e-7, link = "logit", emmax = 50) {
  n <- nrow(fit$X)
  nb <- length(fit$bnm)
  nbeta <- length(fit$betanm)
  censoring_rates <- numeric(nboot)
  percentage_plateau <- numeric(nboot)
  max_event_times <- numeric(nboot)
  
  tempdata <- cbind(fit$Time, fit$Status, fit$X, fit$Z)
  data1 <- subset(tempdata, fit$Status == 1)
  data0 <- subset(tempdata, fit$Status == 0)
  n1 <- nrow(data1)
  n0 <- nrow(data0)
  
  b_boot <- matrix(0, nboot, nb)
  beta_boot <- matrix(0, nboot, nbeta)
  
  i <- 1
  while (i <= nboot) {
    #print(paste0("Round:", i))
    id1 <- sample(1:n1, n1, replace = TRUE)
    id0 <- sample(1:n0, n0, replace = TRUE)
    bootdata <- rbind(data1[id1,], data0[id0,])
    
    bootZ <- bootdata[, fit$bnm]
    bootX <- as.matrix(cbind(rep(1, n), bootdata[, fit$betanm]))
    
    bootfit <- em(bootdata[, 1], bootdata[, 2], bootX, bootZ, NULL, fit$b, fit$beta, link, emmax, eps)
    
    if (bootfit$tau < eps) {
      b_boot[i, ] <- bootfit$b
      beta_boot[i, ] <- bootfit$latencyfit
      Delta_star <- bootdata[, 2]
      Y_star <- bootdata[, 1]
      censoring_rates[i] <- 1 - mean(Delta_star)
      
      last_event_time <- max(Y_star[Delta_star == 1])
      max_event_times[i] <- last_event_time
      censored_after <- sum(Y_star[Delta_star == 0] > last_event_time)
      percentage_plateau[i] <- censored_after / n
      
      i <- i + 1
    }
  }
  
  b_var <- apply(b_boot, 2, var)
  beta_var <- apply(beta_boot, 2, var)
  b_sd <- sqrt(b_var)
  beta_sd <- sqrt(beta_var)
  
  list(b_boot = b_boot, beta_boot = beta_boot, b_sd = b_sd, beta_sd = beta_sd, censoring_rates = censoring_rates, percentage_plateau = percentage_plateau, max_event_times = max_event_times, error = bootfit$error)
}


## ------------------------------------------------------------------------------------------------------------------------------------
# Conditional bootstrap algorithm 1
conditional_bootstrap_1 <- function(fit, nboot = 500, eps = 1e-7, link = "logit", emmax = 50) {
  n <- nrow(fit$X)
  nb <- length(fit$bnm)
  nbeta <- length(fit$betanm)
  censoring_rates <- numeric(nboot)
  percentage_plateau <- numeric(nboot)
  max_event_times <- numeric(nboot)
  
  # Return baseline cumulative hazard (for T*)
  death_times_unique <- sort(unique(fit$Time[fit$Status == 1]))
  cumhaz_values <- fit$cumhaz
  
  # Estimate Kaplan-Meier survival for censoring (for C*)
  censor_status <- 1 - fit$Status
  km_censoring <- survfit(Surv(fit$Time, censor_status) ~ 1)
  censor_times <- km_censoring$time
  censor_surv <- km_censoring$surv
  
  # Remove duplicated censoring survival values
  unique_idx_cens <- !duplicated(1 - censor_surv)
  censor_times_unique <- censor_times[unique_idx_cens]
  censor_surv_unique <- censor_surv[unique_idx_cens]
  
  b_boot <- matrix(0, nboot, nb)
  beta_boot <- matrix(0, nboot, nbeta)
  
  i <- 1
  while (i <= nboot) {
    # New bootstrap dataset
    T_star <- numeric(n)
    C_star <- numeric(n)
    Y_star <- numeric(n)
    Delta_star <- numeric(n)
    
    for (j in 1:n) {
      # Calculate uncure probability
      pi <- exp(crossprod(fit$b, fit$Z[j,])) / (1 + exp(crossprod(fit$b, fit$Z[j,])))
      
      if (runif(1) <= pi) {
        # Not cured: simulate event time
        U <- runif(1)
        target_hazard <- - (1 / exp(as.numeric(crossprod(fit$beta, fit$X[j, -1])))) * log(U)
        idx <- which(cumhaz_values >= target_hazard)[1]
        T_star[j] <- if (!is.na(idx)) death_times_unique[idx] else max(death_times_unique)
        #}
      } else {
        # Cured
        T_star[j] <- Inf
      }
      
      # Simulate censoring time
      U_cens <- runif(1)
      idx_cens <- which(censor_surv_unique <= U_cens)[1]
      C_star[j] <- censor_times_unique[idx_cens]
      
      # Observed outcomes
      Y_star[j] <- min(T_star[j], C_star[j])
      Delta_star[j] <- as.integer(T_star[j] <= C_star[j])
    }
    
    censoring_rates[i] <- 1 - mean(Delta_star)
    
    km_boot <- survfit(Surv(Y_star, Delta_star) ~ 1)
    last_event_time_boot <- max(Y_star[Delta_star == 1])
    max_event_times[i] <- last_event_time_boot
    last_event_time_true <- max(fit$Time[fit$Status == 1])
    censored_times_boot <- Y_star[Delta_star == 0]
    plateau_censored_boot <- censored_times_boot[censored_times_boot > last_event_time_boot]
    percentage_plateau[i] <- length(plateau_censored_boot)/n
    
    # Create bootstrap sample
    bootdata <- cbind(Y_star, Delta_star, fit$X, fit$Z)
    bootZ <- bootdata[, fit$bnm]
    bootX <- as.matrix(cbind(rep(1, n), bootdata[, fit$betanm]))
    
    bootfit <- em(Y_star, Delta_star, bootX, bootZ, NULL, fit$b, fit$beta, link, emmax, eps)
    
    if (bootfit$tau < eps) {
      b_boot[i, ] <- bootfit$b
      beta_boot[i, ] <- bootfit$latencyfit
      i <- i + 1
    }
  }
  
  b_var <- apply(b_boot, 2, var)
  beta_var <- apply(beta_boot, 2, var)
  b_sd <- sqrt(b_var)
  beta_sd <- sqrt(beta_var)
  
  list(b_boot = b_boot, beta_boot = beta_boot, b_sd = b_sd, beta_sd = beta_sd, censoring_rates = censoring_rates, percentage_plateau = percentage_plateau, max_event_times = max_event_times, error = bootfit$error)
}

## ------------------------------------------------------------------------------------------------------------------------------------
# Conditional bootstrap algorithm 2
conditional_bootstrap_2 <- function(fit, nboot = 500, eps = 1e-7, link = "logit", emmax = 50) {
  n <- nrow(fit$X)
  nb <- length(fit$bnm)
  nbeta <- length(fit$betanm)
  censoring_rates <- numeric(nboot)
  percentage_plateau <- numeric(nboot)
  max_event_times <- numeric(nboot)
  
  # Return baseline cumulative hazard (for T*)
  death_times_unique <- sort(unique(fit$Time[fit$Status == 1]))
  cumhaz_values <- fit$cumhaz
  
  # Estimate Kaplan-Meier survival for censoring (for C*)
  censor_status <- 1 - fit$Status
  km_censoring <- survfit(Surv(fit$Time, censor_status) ~ 1)
  censor_times <- km_censoring$time
  censor_surv <- km_censoring$surv
  
  # Remove duplicated censoring survival values
  unique_idx_cens <- !duplicated(1 - censor_surv)
  censor_times_unique <- censor_times[unique_idx_cens]
  censor_surv_unique <- censor_surv[unique_idx_cens]
  
  b_boot <- matrix(0, nboot, nb)
  beta_boot <- matrix(0, nboot, nbeta)
  
  i <- 1
  while (i <= nboot) {
    T_star <- numeric(n)
    C_star <- numeric(n)
    Y_star <- numeric(n)
    Delta_star <- numeric(n)
    
    for (j in 1:n) {
      # Calculate uncure probability
      pi <- exp(crossprod(fit$b, fit$Z[j,])) / (1 + exp(crossprod(fit$b, fit$Z[j,])))
      
      if (runif(1) <= pi) {
        # Not cured: simulate event time
        U <- runif(1)
        target_hazard <- - (1 / exp(as.numeric(crossprod(fit$beta, fit$X[j, -1])))) * log(U)
        idx <- which(cumhaz_values >= target_hazard)[1]
        T_star[j] <- if (!is.na(idx)) death_times_unique[idx] else max(death_times_unique)
      } else {
        # Cured
        T_star[j] <- Inf
      }
      
      # Resample censoring time
      if (fit$Status[j] == 0) {
        # Exact observed censoring
        C_star[j] <- fit$Time[j]
      } else {
        # Conditional censoring distribution: resample for Delta = 1
        idx_yi <- max(which(censor_times_unique <= fit$Time[j]))
        G_Yi <- censor_surv_unique[idx_yi]
        
        U_cens <- runif(1)
        G_target <- G_Yi - G_Yi * U_cens
        
        idx_cens <- which(censor_surv_unique <= G_target)[1]
        C_star[j] <- censor_times_unique[idx_cens]
      }
      
      # Observed outcomes
      Y_star[j] <- min(T_star[j], C_star[j])
      Delta_star[j] <- as.integer(T_star[j] <= C_star[j])
    }
    censoring_rates[i] <- 1-mean(Delta_star)
    
    km_boot <- survfit(Surv(Y_star, Delta_star) ~ 1)
    last_event_time_boot <- max(Y_star[Delta_star == 1])
    max_event_times[i] <- last_event_time_boot
    last_event_time_true <- max(fit$Time[fit$Status == 1])
    censored_times_boot <- Y_star[Delta_star == 0]
    plateau_censored_boot <- censored_times_boot[censored_times_boot > last_event_time_boot]
    percentage_plateau[i] <- length(plateau_censored_boot)/n
    
    bootdata <- cbind(Y_star, Delta_star, fit$X, fit$Z)
    bootZ <- bootdata[, fit$bnm]
    bootX <- as.matrix(cbind(rep(1, n), bootdata[, fit$betanm]))
    
    bootfit <- em(Y_star, Delta_star, bootX, bootZ, NULL, fit$b, fit$beta, link, emmax, eps)
    
    if (bootfit$tau < eps) {
      b_boot[i, ] <- bootfit$b
      beta_boot[i, ] <- bootfit$latencyfit
      i <- i + 1
    }
  }
  
  b_var <- apply(b_boot, 2, var)
  beta_var <- apply(beta_boot, 2, var)
  b_sd <- sqrt(b_var)
  beta_sd <- sqrt(beta_var)
  
  list(b_boot = b_boot, beta_boot = beta_boot, b_sd = b_sd, beta_sd = beta_sd, censoring_rates = censoring_rates, percentage_plateau = percentage_plateau, max_event_times = max_event_times, error = bootfit$error)
} 


## ------------------------------------------------------------------------------------------------------------------------------------
# EM algorithm (PH model only)
em <- function(Time, Status, X, Z, offsetvar, b, beta, link, emmax, eps) {
  error <- FALSE
  w <- Status
  n <- length(Status)
  
  s <- smsurv(Time, Status, X, beta, w)$survival
  HHazard <- smsurv(Time, Status, X, beta, w)$HHazard
  cumhaz <- smsurv(Time, Status, X, beta, w)$cumhaz
  
  convergence <- 1000
  i <- 1
  
  while (convergence > eps & i < emmax) {
    uncureprob <- matrix(exp(Z %*% b) / (1 + exp(Z %*% b)), ncol = 1)
    survival <- drop(s^(exp(beta %*% t(X[,-1]))))
    
    if (any(is.na(uncureprob))) {
      error <- TRUE
      na_idx <- which(is.na(uncureprob))
      break
    }
    if (any(is.na(survival))) {
      error <- TRUE
      print("survival has NA")
      break
    }
    
    w <- Status + (1 - Status) * (uncureprob * survival) / ((1 - uncureprob) + uncureprob * survival)
    
    last_event_time <- max(Time[Status == 1])
    plateau_indices <- which(Status == 0 & Time > last_event_time)
    w[plateau_indices] <- 0
    
    logistfit <- glm(w ~ Z[,-1], family = quasibinomial(link = link))
    update_cureb <- logistfit$coef
    
    if (!is.null(offsetvar)) {
      logistfit <- glm(w ~ Z[,-1] + offset(offsetvar), family = quasibinomial(link = link))
      update_cureb <- logistfit$coef
    }
    
    update_beta <- coxph(Surv(Time, Status) ~ X[,-1] + offset(log(w)), subset = w != 0, method = "breslow")$coef
    
    if (!is.null(offsetvar)) {
      update_beta <- coxph(Surv(Time, Status) ~ X[,-1] + offset(offsetvar + log(w)), subset = w != 0, method = "breslow")$coef
    }
    
    update_s <- smsurv(Time, Status, X, beta, w)$survival
    
    convergence <- sum(c(update_cureb - b, update_beta - beta)^2) + sum((s - update_s)^2)
    
    b <- update_cureb
    beta <- update_beta
    s <- update_s
    
    i <- i + 1
  }
  
  list(logistfit = logistfit, b = b, latencyfit = beta, Survival = s, HHazard = HHazard, cumhaz = cumhaz, Uncureprob = uncureprob, tau = convergence, error = error)
}
  
## ------------------------------------------------------------------------------------------------------------------------------------
# Main SMCURE function - ONLY fitting, no bootstrap
smcure <- function(formula, cureform, offset = NULL, data, na.action = na.omit, 
                   link = "logit", emmax = 50, eps = 1e-7) {
  call <- match.call()
  cat("Program is running..be patient...")
  
  ## Prepare data
  data <- na.action(data)
  n <- nrow(data)
  
  mf <- model.frame(formula, data)
  cvars <- all.vars(cureform)
  Z <- as.matrix(cbind(rep(1, n), data[, cvars]))
  colnames(Z) <- c("(Intercept)", cvars)
  
  if (!is.null(offset)) {
    offsetvar <- data[, all.vars(offset)]
  } else {
    offsetvar <- NULL
  }
  
  Y <- model.extract(mf, "response")
  X <- model.matrix(attr(mf, "terms"), mf)
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  
  Time <- Y[, 1]
  Status <- Y[, 2]
  
  bnm <- colnames(Z)
  nb <- ncol(Z)
  betanm <- colnames(X)[-1]
  nbeta <- ncol(X) - 1
  
  ## Initial values
  w <- Status
  b <- glm(w ~ Z[,-1], family = quasibinomial(link = link))$coef
  beta <- coxph(Surv(Time, Status) ~ X[,-1] + offset(log(w)), subset = w != 0, method = "breslow")$coef
  
  ## EM algorithm
  emfit <- em(Time, Status, X, Z, offsetvar, b, beta, link, emmax, eps)
  
  fit <- list()
  class(fit) <- "smcure"
  
  fit$logistfit <- emfit$logistfit
  fit$b <- emfit$b
  fit$beta <- emfit$latencyfit
  fit$s <- emfit$Survival
  fit$HHazard <- emfit$HHazard
  fit$cumhaz <- emfit$cumhaz
  fit$Time <- Time
  fit$Status <- Status
  fit$X <- X
  fit$Z <- Z
  fit$bnm <- bnm
  fit$betanm <- betanm
  fit$call <- call
  fit$error <- emfit$error
  
  cat(" done.\n")
  
  return(fit)
}

## ------------------------------------------------------------------------------------------------------------------------------------
# Survival function (PH model only)
smsurv <- function(Time, Status, X, beta, w) {
  death_point <- sort(unique(subset(Time, Status == 1)))
  
  coxexp <- exp(beta %*% t(X[,-1]))
  
  lambda <- numeric()
  event <- numeric()
  
  for (i in 1:length(death_point)) {
    event[i] <- sum(Status * as.numeric(Time == death_point[i]))
    temp <- sum(as.numeric(Time >= death_point[i]) * w * drop(coxexp))
    temp1 <- event[i]
    lambda[i] <- temp1 / temp
  }
  
  cumhaz_values <- cumsum(lambda)
  
  HHazard <- numeric()
  
  if (anyNA(Time)) stop("Time vector contains NA values")
  if (anyNA(Status)) stop("Status vector contains NA values")
  
  for (i in 1:length(Time)) {
    HHazard[i] <- sum(as.numeric(Time[i] >= death_point) * lambda)
    if (Time[i] > max(death_point)) HHazard[i] <- Inf
    if (Time[i] < min(death_point)) HHazard[i] <- 0
  }
  
  survival <- exp(-HHazard)
  
  list(survival = survival, HHazard = HHazard, cumhaz = cumhaz_values)
}

## ------------------------------------------------------------------------------------------------------------------------------------
# Setup
N <- 500
beta <- c(1, -0.1, 0.8)
gamma <- c(-0.8, 1.3, 1.5, -0.2)
alpha <- 0.05
nboot <- 200
methods <- c("naive", "naive stratified", "conditional1", "conditional2")

# Initialize storage vector
init_vec <- function() rep(NA, N)

# Initialize an empty list to store results per method
results_list <- list()

# Initialize method-specific structures before looping
for (method in methods) {
  results_list[[method]] <- list()
  
  for (ci_type in c("normal", "standard", "percentile")) {
    for (j in 1:4) {
      results_list[[method]][[paste0("coverage_", ci_type, "_b", j)]] <- init_vec()
      results_list[[method]][[paste0("length_", ci_type, "_b", j)]] <- init_vec()
    }
    for (j in 1:3) {
      results_list[[method]][[paste0("coverage_", ci_type, "_beta", j)]] <- init_vec()
      results_list[[method]][[paste0("length_", ci_type, "_beta", j)]] <- init_vec()
    }
  }
  
  results_list[[method]]$censoring_rates_boot <- init_vec()
  results_list[[method]]$percentage_plateau_boot <- init_vec()
  results_list[[method]]$max_event_times_boot <- init_vec()
}

censoring_full_vec <- init_vec()
plateau_full_vec <- init_vec()
max_event_time_full_vec <- init_vec()

# Outer loop over repetitions
for (i in 1:N) {
  repeat {
    simulated_data <- simulate_data(n = 500)
    
    fit <- smcure(Surv(simulated_data$time, simulated_data$censored) ~ X.1 + X.2 + X.3,
                  cureform = ~ Z.2 + Z.3 + Z.4, data = simulated_data)
    
    if (fit$error == TRUE) {
      print("An error occured, resampling...")
      next
    }
    
    for (method in methods) {
      print(paste0("progress: ", round(((i - 1) * length(methods) + match(method, methods)) / (N * length(methods)), 4) * 100))
      
      any_bootstrap_failed <- FALSE
      results <- do_bootstrap(fit, method = method, nboot = nboot)
      if (!is.null(results$error) && results$error == TRUE) {
        print(paste0("Bootstrap error with method ", method, ". Resampling dataset..."))
        any_bootstrap_failed <- TRUE
        break  # Exit this inner method loop
      }
      
      fit$b_sd <- results$b_sd
      fit$beta_sd <- results$beta_sd
      
      z_crit <- qnorm(1 - alpha / 2)
      
      # Normal CI
      ci_normal_lower_b <- fit$b - z_crit * fit$b_sd
      ci_normal_upper_b <- fit$b + z_crit * fit$b_sd
      
      for (j in 1:4) {
        results_list[[method]][[paste0("coverage_normal_b", j)]][i] <- as.numeric(gamma[j] >= ci_normal_lower_b[j] & gamma[j] <= ci_normal_upper_b[j])
        results_list[[method]][[paste0("length_normal_b", j)]][i] <- 2 * z_crit * fit$b_sd[j]
      }
      
      ci_normal_lower_beta <- fit$beta - z_crit * fit$beta_sd
      ci_normal_upper_beta <- fit$beta + z_crit * fit$beta_sd
      
      for (j in 1:3) {
        results_list[[method]][[paste0("coverage_normal_beta", j)]][i] <- as.numeric(beta[j] >= ci_normal_lower_beta[j] & beta[j] <= ci_normal_upper_beta[j])
        results_list[[method]][[paste0("length_normal_beta", j)]][i] <- 2 * z_crit * fit$beta_sd[j]
      }
      
      # Standard CI
      n_b <- nrow(results$b_boot)
      n_beta <- nrow(results$beta_boot)
      
      b_diff <- sweep(results$b_boot, 2, fit$b)
      beta_diff <- sweep(results$beta_boot, 2, fit$beta)
      
      ci_standard_lower_b <- fit$b - (n_b^(-1/2)) * apply(sqrt(n_b) * b_diff, 2, quantile, probs = 1 - alpha / 2)
      ci_standard_upper_b <- fit$b - (n_b^(-1/2)) * apply(sqrt(n_b) * b_diff, 2, quantile, probs = alpha / 2)
      
      for (j in 1:4) {
        results_list[[method]][[paste0("coverage_standard_b", j)]][i] <- as.numeric(gamma[j] >= ci_standard_lower_b[j] & gamma[j] <= ci_standard_upper_b[j])
        results_list[[method]][[paste0("length_standard_b", j)]][i] <- ci_standard_upper_b[j] - ci_standard_lower_b[j]
      }
      
      ci_standard_lower_beta <- fit$beta - (n_beta^(-1/2)) * apply(sqrt(n_beta) * beta_diff, 2, quantile, probs = 1 - alpha / 2)
      ci_standard_upper_beta <- fit$beta - (n_beta^(-1/2)) * apply(sqrt(n_beta) * beta_diff, 2, quantile, probs = alpha / 2)
      
      for (j in 1:3) {
        results_list[[method]][[paste0("coverage_standard_beta", j)]][i] <- as.numeric(beta[j] >= ci_standard_lower_beta[j] & beta[j] <= ci_standard_upper_beta[j])
        results_list[[method]][[paste0("length_standard_beta", j)]][i] <- ci_standard_upper_beta[j] - ci_standard_lower_beta[j]
      }
      
      # Percentile CI
      ci_percentile_lower_b <- apply(results$b_boot, 2, quantile, probs = alpha / 2)
      ci_percentile_upper_b <- apply(results$b_boot, 2, quantile, probs = 1 - alpha / 2)
      
      for (j in 1:4) {
        results_list[[method]][[paste0("coverage_percentile_b", j)]][i] <- as.numeric(gamma[j] >= ci_percentile_lower_b[j] & gamma[j] <= ci_percentile_upper_b[j])
        results_list[[method]][[paste0("length_percentile_b", j)]][i] <- ci_percentile_upper_b[j] - ci_percentile_lower_b[j]
      }
      
      ci_percentile_lower_beta <- apply(results$beta_boot, 2, quantile, probs = alpha / 2)
      ci_percentile_upper_beta <- apply(results$beta_boot, 2, quantile, probs = 1 - alpha / 2)
      
      for (j in 1:3) {
        results_list[[method]][[paste0("coverage_percentile_beta", j)]][i] <- as.numeric(beta[j] >= ci_percentile_lower_beta[j] & beta[j] <= ci_percentile_upper_beta[j])
        results_list[[method]][[paste0("length_percentile_beta", j)]][i] <- ci_percentile_upper_beta[j] - ci_percentile_lower_beta[j]
      }
      
      results_list[[method]]$censoring_rates_boot[i] <- mean(results$censoring_rates)
      results_list[[method]]$percentage_plateau_boot[i] <- mean(results$percentage_plateau)
      results_list[[method]]$max_event_times_boot[i] <- mean(results$max_event_times)
    }
    
    if (any_bootstrap_failed) {
      next  # Retry outer repeat (i-th iteration)
    }
    
    break
    
  }
  
  last_event_time <- max(simulated_data$time[simulated_data$censored == 1])
  plateau_censored <- simulated_data$time[simulated_data$censored == 0 & simulated_data$time > last_event_time]
  plateau_full_vec[i] <- length(plateau_censored) / nrow(simulated_data)
  censoring_full_vec[i] <- 1 - mean(simulated_data$censored)
  max_event_time_full_vec[i] <- last_event_time
  
}

# Summarize coverage and lengths
param_summary <- do.call(rbind, lapply(names(results_list), function(method) {
  res <- results_list[[method]]
  data.frame(
    Method = rep(method, 7),
    Parameter = c(paste0("b", 1:4), paste0("beta", 1:3)),
    True_Value = c(gamma, beta),
    Coverage_Standard_CI = c(sapply(1:4, function(j) mean(res[[paste0("coverage_standard_b", j)]])), sapply(1:3, function(j) mean(res[[paste0("coverage_standard_beta", j)]]))),
    Coverage_Normal_CI = c(sapply(1:4, function(j) mean(res[[paste0("coverage_normal_b", j)]])), sapply(1:3, function(j) mean(res[[paste0("coverage_normal_beta", j)]]))),
    Coverage_Percentile_CI = c(sapply(1:4, function(j) mean(res[[paste0("coverage_percentile_b", j)]])), sapply(1:3, function(j) mean(res[[paste0("coverage_percentile_beta", j)]]))),
    Length_Standard_CI = c(sapply(1:4, function(j) mean(res[[paste0("length_standard_b", j)]])), sapply(1:3, function(j) mean(res[[paste0("length_standard_beta", j)]]))),
    Length_Normal_CI = c(sapply(1:4, function(j) mean(res[[paste0("length_normal_b", j)]])), sapply(1:3, function(j) mean(res[[paste0("length_normal_beta", j)]]))),
    Length_Percentile_CI = c(sapply(1:4, function(j) mean(res[[paste0("length_percentile_b", j)]])), sapply(1:3, function(j) mean(res[[paste0("length_percentile_beta", j)]])))
  )
}))

param_summary <- param_summary[order(param_summary$Parameter, param_summary$Method), ]
print(param_summary)

# Summarize censoring and plateau
censor_plateau_summary <- data.frame(
  Source = c("Full dataset", paste0(names(results_list), " bootstrap average")),
  Censoring_Percentage = c(
    mean(censoring_full_vec),
    sapply(results_list, function(r) mean(r$censoring_rates_boot))
  ),
  Plateau_Percentage = c(
    mean(plateau_full_vec),
    sapply(results_list, function(r) mean(r$percentage_plateau_boot))
  ),
  Max_event_time = c(
    mean(max_event_time_full_vec),
    sapply(results_list, function(r) mean(r$max_event_times_boot))
  )
)

print(censor_plateau_summary)

