# ============================================================
# Main program: common libraries, utilities, fitters, selectors, and data generators
# ============================================================

# Optimized Simulation Studies for PMHD-MCP AEPD
# CORRECTED VERSION 鈥?Key changes from original:
#   1. Added AEPD-MLE competitor (EM for AEPD mixtures)
#   2. Fixed Hellinger distance fairness: added ISE metric +
#      report K selected by each method so readers can see
#      when GM gains come from inflated K
#   3. Elbow rule thresholds aligned with paper (0.70, 0.20, 2.0)
#   4. Scenario 3: skewed-t alpha_skew set to 5 (was 0 = symmetric)
#   5. Scenario 4: added GM-MLE and AEPD-MLE competitors
#   6. Scenario 2: fixed H_gmmle_K0 extraction bug
#   7. Scenario 5: added QQ-plot generation
#   8. Alpha search bounds documented
# ============================================================

library(MASS)
library(nloptr)
library(mclust)
library(clue)
library(foreach)
library(doParallel)

# --- Parallel backend helpers ---
start_parallel_backend <- function(parallel, ncores) {
  if (!parallel) {
    foreach::registerDoSEQ()
    return(NULL)
  }
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  cl
}

stop_parallel_backend <- function(cl) {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
  }
  foreach::registerDoSEQ()
  invisible(gc())
}

# Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- Core Utilities ---
trapz_weights <- function(x) {
  n <- length(x); dx <- diff(x); w <- numeric(n)
  w[1] <- dx[1] / 2; w[n] <- dx[n - 1] / 2
  if (n > 2) w[2:(n - 1)] <- (dx[1:(n - 2)] + dx[2:(n - 1)]) / 2
  w
}

normalize_grid_density <- function(y, w, eps = 1e-15) {
  y <- pmax(y, 0); z <- sum(y * w)
  if (!is.finite(z) || z <= eps) return(rep(1/length(y), length(y)))
  y / z
}

H2_grid <- function(p, q, w) {
  p <- pmax(p, 0); q <- pmax(q, 0)
  val <- 1 - sum(sqrt(p * q) * w)
  2 * pmax(val, 0)
}

# --- ISE (Integrated Squared Error) ---
# Model-agnostic metric: ISE = integral (f_hat - f_0)^2 dx
# Unlike Hellinger, ISE does not structurally favor the estimation
# criterion of any particular method.
ISE_grid <- function(p, q, w) {
  sum((p - q)^2 * w)
}

mcp_value <- function(pi, lambda, gamma)
  sum(ifelse(pi <= gamma * lambda,
             lambda * pi - pi^2 / (2 * gamma),
             gamma * lambda^2 / 2))

mcp_prime <- function(pi, lambda, gamma) pmax(lambda - pi / gamma, 0)

# --- AEPD Density ---
daepd <- function(x, mu, sigma, alpha, tau) {
  sigma <- pmax(sigma, 1e-10); alpha <- pmax(alpha, 1e-10)
  tau   <- pmin(pmax(tau, 1e-8), 1 - 1e-8)
  C <- alpha * tau * (1 - tau) / (gamma(1 / alpha) * sigma)
  z <- abs(x - mu)^alpha / sigma^alpha
  d <- C * exp(-z * ifelse(x < mu, (1 - tau)^alpha, tau^alpha))
  pmax(d, 1e-300)
}

# --- Initialization ---
init_params_kmeans <- function(x, K, alpha_init = 2.0) {
  km  <- kmeans(matrix(x, ncol = 1), centers = K, nstart = 20)
  ord <- order(km$centers[, 1]); mu <- sort(km$centers[, 1])
  cl_new <- match(km$cluster, ord)
  pi <- tabulate(cl_new, nbins = K) / length(x)
  global_sd <- sd(x)
  sigma <- sapply(seq_len(K), function(k) {
    xs <- x[cl_new == k]
    s  <- suppressWarnings(sd(xs))
    if (!is.finite(s) || s < 0.1 * global_sd) s <- global_sd / 3
    s
  })
  list(pi = pi / sum(pi), mu = mu,
       sigma = pmax(sigma, 0.1 * global_sd),
       alpha = rep(alpha_init, K), tau = rep(0.5, K))
}

# --- PMHD-MCP AEPD Core Fitter ---
pmhd_mcp_fit_one <- function(x, K_init, lambda, gamma = 3,
                             grid, wq, ghat,
                             tol = 1e-5, max_iter = 100,
                             delta_inner = 1e-8,
                             delta_final = 1e-3, init = NULL,
                             maxeval = 250,
                             alpha_lower = 1.0, alpha_upper = 3.0,
                             tau_lower = 0.01, tau_upper = 0.99) {
  sigma_lower <- 1e-3 * sd(x); sigma_upper <- 3 * diff(range(x))
  log_sl <- log(sigma_lower); log_su <- log(sigma_upper)
  log_al <- log(alpha_lower); log_au <- log(alpha_upper)

  if (is.null(init)) init <- init_params_kmeans(x, K_init)
  pi_k <- pmax(init$pi, delta_inner); pi_k <- pi_k / sum(pi_k)
  mu_k <- init$mu; sigma_k <- init$sigma
  alpha_k <- init$alpha; tau_k <- init$tau
  obj_old <- Inf

  for (iter in seq_len(max_iter)) {
    f_old <- matrix(0, length(grid), K_init)
    for (k in seq_len(K_init))
      f_old[, k] <- daepd(grid, mu_k[k], sigma_k[k], alpha_k[k], tau_k[k])
    f_mix_old <- pmax(drop(f_old %*% pi_k), 1e-300)
    R_old <- sweep(sweep(f_old, 2, pi_k, "*"), 1, f_mix_old, "/")
    c_k <- colSums(ghat * R_old * wq)

    for (k in which(pi_k >= delta_inner)) {
      gtilde_k <- (ghat * R_old[, k]) / pmax(c_k[k], 1e-12)
      obj_theta <- function(par) {
        mu <- par[1]; sig <- exp(par[2])
        alp <- exp(par[3]); t <- plogis(par[4])
        1 - sum(sqrt(gtilde_k * daepd(grid, mu, sig, alp, t)) * wq)
      }
      res <- tryCatch(
        nloptr(c(mu_k[k], log(sigma_k[k]),
                 log(alpha_k[k]), qlogis(tau_k[k])),
               obj_theta,
               lb = c(-Inf, log_sl, log_al, qlogis(tau_lower)),
               ub = c( Inf, log_su, log_au, qlogis(tau_upper)),
               opts = list(algorithm = "NLOPT_LN_BOBYQA",
                           xtol_rel = 1e-6, maxeval = maxeval)),
        error = function(e) NULL)
      if (!is.null(res)) {
        mu_k[k]    <- res$solution[1]
        sigma_k[k] <- exp(res$solution[2])
        alpha_k[k] <- exp(res$solution[3])
        tau_k[k]   <- plogis(res$solution[4])
      }
    }

    f_new <- matrix(0, length(grid), K_init)
    for (k in seq_len(K_init))
      f_new[, k] <- daepd(grid, mu_k[k], sigma_k[k], alpha_k[k], tau_k[k])
    B_k <- sapply(seq_len(K_init), function(k)
      2 * sum(sqrt(pmax(ghat * R_old[, k] * f_new[, k], 0)) * wq))
    xi_k <- mcp_prime(pi_k, lambda, gamma)

    g_nu <- function(nu) sum(pmax(B_k / (2 * (xi_k + nu)), 0)^2) - 1
    lower_nu <- -min(xi_k) + 1e-10
    upper_nu <- max(B_k) + 1

    if (all(xi_k <= 0)) {
      pi_k <- B_k^2
    } else if (g_nu(lower_nu) * g_nu(upper_nu) < 0) {
      nu_star <- uniroot(g_nu, c(lower_nu, upper_nu), tol = 1e-10)$root
      pi_k <- pmax(B_k / (2 * (xi_k + nu_star)), 0)^2
    } else {
      pi_k <- rep(0, K_init); pi_k[which.max(B_k)] <- 1
    }
    if (sum(pi_k) <= 0) pi_k <- rep(1 / K_init, K_init)
    else pi_k <- pi_k / sum(pi_k)

    f_mix_new <- normalize_grid_density(pmax(drop(f_new %*% pi_k), 1e-300), wq)
    H2_now <- H2_grid(ghat, f_mix_new, wq)
    obj_now <- H2_now + mcp_value(pi_k, lambda, gamma)
    if (abs(obj_old - obj_now) < tol) break
    obj_old <- obj_now
  }

  keep <- which(pi_k >= delta_final)
  if (length(keep) == 0L) {
    keep <- which.max(pi_k)
  }
  list(pi = pi_k[keep] / sum(pi_k[keep]),
       mu = mu_k[keep], sigma = sigma_k[keep],
       alpha = alpha_k[keep], tau = tau_k[keep],
       Khat = length(keep), H2 = H2_now)
}


# ============================================================
# AEPD-MLE competitor (NEW 鈥?was missing from original code)
# EM algorithm for AEPD mixture + BIC order selection
# ============================================================

loglik_aepd_mix <- function(x, pi_k, mu_k, sigma_k, alpha_k, tau_k) {
  K <- length(pi_k); n <- length(x)
  lik_mat <- matrix(0, n, K)
  for (k in seq_len(K))
    lik_mat[, k] <- pi_k[k] * daepd(x, mu_k[k], sigma_k[k],
                                      alpha_k[k], tau_k[k])
  sum(log(pmax(rowSums(lik_mat), 1e-300)))
}

em_aepd_mix <- function(x, K, max_iter = 200, tol = 1e-6,
                        alpha_lower = 0.3, alpha_upper = 5.0,
                        tau_lower = 0.01, tau_upper = 0.99,
                        init = NULL, maxeval = 200) {
  n <- length(x)
  sigma_lower <- 1e-3 * sd(x)
  sigma_upper <- 3 * diff(range(x))

  if (is.null(init)) init <- init_params_kmeans(x, K)
  pi_k <- init$pi; mu_k <- init$mu; sigma_k <- init$sigma
  alpha_k <- init$alpha; tau_k <- init$tau
  ll_old <- -Inf

  for (iter in seq_len(max_iter)) {
    lik_mat <- matrix(0, n, K)
    for (k in seq_len(K))
      lik_mat[, k] <- pi_k[k] * daepd(x, mu_k[k], sigma_k[k],
                                        alpha_k[k], tau_k[k])
    row_sum <- pmax(rowSums(lik_mat), 1e-300)
    gamma_mat <- lik_mat / row_sum

    n_k <- colSums(gamma_mat)
    pi_k <- pmax(n_k / n, 1e-8)
    pi_k <- pi_k / sum(pi_k)

    for (k in seq_len(K)) {
      wk <- gamma_mat[, k]
      obj_k <- function(par) {
        mu <- par[1]; sig <- exp(par[2])
        alp <- exp(par[3]); tau <- plogis(par[4])
        -sum(wk * log(pmax(daepd(x, mu, sig, alp, tau), 1e-300)))
      }
      res <- tryCatch(
        nloptr(c(mu_k[k], log(sigma_k[k]),
                 log(alpha_k[k]), qlogis(tau_k[k])),
               obj_k,
               lb = c(-Inf, log(sigma_lower),
                      log(alpha_lower), qlogis(tau_lower)),
               ub = c( Inf, log(sigma_upper),
                      log(alpha_upper), qlogis(tau_upper)),
               opts = list(algorithm = "NLOPT_LN_BOBYQA",
                           xtol_rel = 1e-6, maxeval = maxeval)),
        error = function(e) NULL)
      if (!is.null(res)) {
        mu_k[k]    <- res$solution[1]
        sigma_k[k] <- exp(res$solution[2])
        alpha_k[k] <- exp(res$solution[3])
        tau_k[k]   <- plogis(res$solution[4])
      }
    }

    ll_new <- loglik_aepd_mix(x, pi_k, mu_k, sigma_k, alpha_k, tau_k)
    if (abs(ll_new - ll_old) < tol * (1 + abs(ll_old))) {
      ll_old <- ll_new
      break
    }
    ll_old <- ll_new
  }

  df <- (K - 1) + K * 4
  bic <- -2 * ll_old + df * log(n)
  list(pi = pi_k, mu = mu_k, sigma = sigma_k,
       alpha = alpha_k, tau = tau_k,
       Khat = K, loglik = ll_old, bic = bic, df = df)
}

init_params_quantile <- function(x, K, alpha_lower = 0.3, alpha_upper = 5.0) {
  centers <- as.numeric(stats::quantile(x, probs = ((seq_len(K) - 0.5) / K),
                                        names = FALSE, type = 8))
  cl <- max.col(-abs(outer(x, centers, "-")), ties.method = "random")
  global_sd <- sd(x)
  pi <- tabulate(cl, nbins = K) / length(x)
  mu <- vapply(seq_len(K), function(k) {
    xs <- x[cl == k]
    if (length(xs) == 0L) centers[k] else mean(xs)
  }, numeric(1))
  sigma <- vapply(seq_len(K), function(k) {
    xs <- x[cl == k]
    s <- suppressWarnings(sd(xs))
    if (!is.finite(s) || s < 0.1 * global_sd) s <- global_sd / 3
    s
  }, numeric(1))
  ord <- order(mu)
  list(pi = pmax(pi[ord], 1e-6) / sum(pmax(pi[ord], 1e-6)),
       mu = mu[ord],
       sigma = pmax(sigma[ord], 0.1 * global_sd),
       alpha = runif(K, max(alpha_lower, 0.5), min(alpha_upper, 3.0)),
       tau = runif(K, 0.2, 0.8))
}

split_aepd_component_init <- function(prev, jitter = 0.35) {
  big <- which.max(prev$pi)
  step <- jitter * prev$sigma[big]
  init <- list(
    pi    = c(prev$pi[-big], prev$pi[big] / 2, prev$pi[big] / 2),
    mu    = c(prev$mu[-big], prev$mu[big] - step, prev$mu[big] + step),
    sigma = c(prev$sigma[-big], prev$sigma[big], prev$sigma[big]),
    alpha = c(prev$alpha[-big], prev$alpha[big], prev$alpha[big]),
    tau   = c(prev$tau[-big], prev$tau[big], prev$tau[big]))
  init$pi <- init$pi / sum(init$pi)
  ord <- order(init$mu)
  lapply(init, function(v) v[ord])
}

fit_aepd_mle <- function(x, K_max, nstart = 10,
                         alpha_lower = 0.3, alpha_upper = 5.0,
                         max_iter = 500, maxeval = 300,
                         criterion = c("bic", "hq", "aic"),
                         penalty_mult = 1.0) {
  criterion <- match.arg(criterion)
  best_fits <- list()  # best fit per K
  
  for (K in 1:K_max) {
    best_ll_K <- -Inf; best_fit_K <- NULL
    
    for (s in seq_len(nstart)) {
      init <- NULL
      if (s == 1) {
        # Default k-means init
        init <- NULL
      } else if (K >= 2 && !is.null(best_fits[[K - 1]]) && s %% 4 == 2) {
        init <- split_aepd_component_init(best_fits[[K - 1]],
                                          jitter = runif(1, 0.25, 0.75))
      } else if (s %% 4 == 3) {
        init <- init_params_quantile(x, K, alpha_lower, alpha_upper)
      } else {
        ini <- init_params_kmeans(x, K)
        ini$mu <- ini$mu + rnorm(K, 0, 0.05 * sd(x))
        ini$alpha <- runif(K, max(alpha_lower, 0.5), min(alpha_upper, 3.0))
        ini$tau <- runif(K, 0.2, 0.8)
        init <- ini
      }
      
      fit <- tryCatch(
        em_aepd_mix(x, K, init = init,
                    alpha_lower = alpha_lower,
                    alpha_upper = alpha_upper,
                    max_iter = max_iter,
                    maxeval = maxeval),
        error = function(e) NULL)
      if (!is.null(fit) && is.finite(fit$loglik) && fit$loglik > best_ll_K) {
        best_ll_K <- fit$loglik; best_fit_K <- fit
      }
    }
    best_fits[[K]] <- best_fit_K
  }
  
  scores <- sapply(best_fits, function(f) {
    if (is.null(f) || !is.finite(f$loglik)) return(Inf)
    n <- length(x)
    penalty <- switch(criterion,
      aic = 2,
      hq = 2 * log(log(n)),
      bic = log(n))
    -2 * f$loglik + penalty_mult * penalty * f$df
  })
  best_K <- which.min(scores)
  if (length(best_K) == 0L || !is.finite(scores[best_K])) return(NULL)
  for (K in seq_along(best_fits)) {
    if (!is.null(best_fits[[K]])) best_fits[[K]]$selection_score <- scores[K]
  }
  selected <- best_fits[[best_K]]
  selected$all_fits <- best_fits
  selected$criterion_table <- do.call(rbind, lapply(seq_along(best_fits), function(K) {
    fit <- best_fits[[K]]
    if (is.null(fit)) {
      return(data.frame(K = K, loglik = NA_real_,
                        df = (K - 1) + 4 * K, score = NA_real_))
    }
    data.frame(K = fit$Khat, loglik = fit$loglik,
               df = fit$df, score = scores[K])
  }))
  selected$aepd_path <- list(
    fits = best_fits,
    scores = scores,
    criterion = criterion,
    penalty_mult = penalty_mult,
    K_grid = seq_along(best_fits),
    loglik = sapply(best_fits, function(f) if (is.null(f)) NA_real_ else f$loglik),
    df = sapply(best_fits, function(f) if (is.null(f)) NA_real_ else f$df)
  )
  selected
}


# --- Path post-processing and H-elbow selector ---
weighted_lm_sse <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & (w > 0)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  if (length(y) <= 1L) return(0)
  fit <- tryCatch(stats::lm.wfit(cbind(1, x), y, w = w),
                  error = function(e) NULL)
  if (is.null(fit)) return(Inf)
  sum(w * fit$residuals^2)
}

compress_lambda_path <- function(lambda_grid, K_mono, H2_vec, H_vec,
                                 K_raw = NULL) {
  K_levels <- sort(unique(K_mono[is.finite(K_mono)]), decreasing = TRUE)
  if (length(K_levels) == 0L) {
    return(data.frame(K = numeric(0), W = integer(0),
                      lam_range = numeric(0),
                      best_idx = integer(0),
                      lambda = numeric(0), H2_star = numeric(0),
                      H_star = numeric(0), logH_star = numeric(0)))
  }
  if (is.null(K_raw)) K_raw <- K_mono
  lambda_pos <- lambda_grid[is.finite(lambda_grid) & lambda_grid > 0]
  lambda_floor <- if (length(lambda_pos) > 0L) min(lambda_pos) / 2 else 1e-12
  lambda_safe <- pmax(lambda_grid, lambda_floor)
  log_lambda <- log(lambda_safe)
  rows <- lapply(K_levels, function(Kv) {
    idx <- which(K_mono == Kv & K_raw == Kv & is.finite(H2_vec))
    if (length(idx) == 0L) return(NULL)
    best_idx <- idx[which.min(H2_vec[idx])]
    log_lams_at_K <- log_lambda[idx]
    data.frame(K = Kv, W = length(idx),
               lam_range = max(log_lams_at_K) - min(log_lams_at_K),
               best_idx = best_idx,
               lambda = lambda_grid[best_idx],
               H2_star = H2_vec[best_idx],
               H_star = H_vec[best_idx],
               logH_star = log(pmax(H_vec[best_idx], 1e-12)))
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L)
    return(data.frame(K = numeric(0), W = integer(0),
                      lam_range = numeric(0),
                      best_idx = integer(0),
                      lambda = numeric(0), H2_star = numeric(0),
                      H_star = numeric(0), logH_star = numeric(0)))
  do.call(rbind, rows)
}

# FIX: Elbow thresholds now match paper Section 3
# Paper says: c_j >= 0.70 * max_k c_k,  imp_j >= 0.20
# Original had: cross_frac=0.85, imp_min=0.30, flat_ratio=2.5
# ============================================================
# K selection based on MCP path analysis
#
# The MCP path K_mono(位) is a decreasing step function.
# The correct K鈧€ is identified by two complementary signals:
#
#   1. LARGEST DROP (螖K): When MCP removes all redundant 
#      components at once, the biggest K jump occurs.
#
#   2. LONGEST PLATEAU: The correct K persists across many 位 
#      values (true components resist the MCP penalty).
#
#   3. BIGGEST H JUMP (for tie-breaking): When K drops below 
#      K鈧€, H* increases sharply because a true component is lost.
#
# Primary rule:
#   score(K) = 螖K_into 脳 W(K)
#   K_opt = argmax score(K), prefer smaller K among ties
#
# When all 螖K = 1 (stepwise path like 5鈫?鈫?鈫?鈫?):
#   score degenerates to W(K). Combined with the H* jump after 
#   K drops, this identifies the "last K worth keeping".
# ============================================================

choose_K_from_path <- function(K_mono, compressed_df, n_obs = NULL,
                               jump_frac = 0.35, jump_min = 1.0) {
  if (is.null(compressed_df) || nrow(compressed_df) == 0)
    return(list(K_opt = NA_integer_, best_idx = NA_integer_,
                reason = "empty_path", scores = NULL))
  
  compressed_df <- compressed_df[is.finite(compressed_df$H_star), , drop = FALSE]
  J <- nrow(compressed_df)
  if (J == 0L)
    return(list(K_opt = NA_integer_, best_idx = NA_integer_,
                reason = "empty_path", scores = NULL))
  if (J == 1L)
    return(list(K_opt = compressed_df$K[1],
                best_idx = compressed_df$best_idx[1],
                reason = "single_K", scores = NULL))
  
  Kvals <- compressed_df$K
  Wvals <- compressed_df$W
  Rvals <- compressed_df$lam_range
  H2vals <- compressed_df$H2_star
  
  # n multiplier (if not provided, use 1 鈥?nH2 and H2 give same ranking)
  nn <- if (!is.null(n_obs)) n_obs else 1
  nH2vals <- nn * H2vals
  
  delta_K <- c(0, Kvals[-J] - Kvals[-1])
  
  # nH2 jump: how much nH2 increases when K drops by 1
  # This is the core signal: losing a TRUE component causes a large,
  # irreducible jump in the objective; losing a REDUNDANT component 
  # causes nearly zero jump.
  #nH2_jump <- c(pmax(nH2vals[-1] - nH2vals[-J], 0), 0)
  nH2_last <- if (Kvals[J] > 1) nH2vals[J] else 0
  nH2_jump <- c(pmax(nH2vals[-1] - nH2vals[-J], 0), nH2_last)
  
  # Normalized lam_range as stability auxiliary
  R_max <- max(Rvals, na.rm = TRUE)
  R_norm <- if (R_max > 1e-10) Rvals / R_max else rep(0, J)
  
  # Weighted score kept as a fallback/diagnostic. The primary rule is
  # the first substantial nH2 jump: once dropping below K causes a clear
  # loss, keep that largest such K rather than rewarding later underfit
  # plateaus.
  nH2_max <- max(nH2_jump, na.rm = TRUE)
  nH2_norm <- if (nH2_max > 1e-10) nH2_jump / nH2_max else rep(0, J)
  
  scores <- 3 * nH2_norm + R_norm
  scores[1] <- -Inf  # exclude K_init overfit plateau
  #scores[J] <- -Inf  # no outgoing jump from the last plateau
  if (Kvals[J] <= 1) scores[J] <- -Inf  # 鍙湪K=1鏃舵墠鎺掗櫎
  
  
  jump_cut <- max(jump_min, jump_frac * nH2_max)
  candidates <- which(nH2_jump >= jump_cut)
  candidates <- candidates[candidates > 1L]
  if (length(candidates) > 0L) {
    pick <- candidates[which.max(nH2_jump[candidates])]
    reason <- "max_substantial_nH2_jump"
  } else {
    finite_scores <- is.finite(scores)
    if (any(finite_scores)) {
      best_score <- max(scores[finite_scores], na.rm = TRUE)
      tied <- which(abs(scores - best_score) < 1e-10)
      pick <- tied[1]  # largest K among ties (first in decreasing order)
      reason <- "nH2_jump_score_fallback"
    } else {
      pick <- J
      reason <- "last_noninitial_plateau_fallback"
    }
  }
  
  list(K_opt = Kvals[pick],
       best_idx = compressed_df$best_idx[pick],
       reason = reason,
       scores = data.frame(K = Kvals, W = Wvals,
         lam_range = round(Rvals, 5),
         delta_K = delta_K,
         nH2 = round(nH2vals, 4),
         nH2_jump = round(nH2_jump, 4),
         jump_cut = round(jump_cut, 4),
         substantial_jump = nH2_jump >= jump_cut,
         R_norm = round(R_norm, 4),
         score = round(scores, 4)))
}

choose_H_elbow <- function(compressed_df,
                           cross_frac = 0.70,
                           imp_min    = 0.20,
                           flat_ratio = 2.0,
                           n_obs = NULL,
                           jump_frac = 0.35,
                           jump_min = 1.0) {
  res <- choose_K_from_path(NULL, compressed_df, n_obs = n_obs,
                            jump_frac = jump_frac, jump_min = jump_min)
  list(best_idx    = res$best_idx,
       K_opt       = res$K_opt,
       split_idx   = NA_integer_,
       fallback    = (res$reason == "single_K" || res$reason == "empty_path"),
       reason      = res$reason,
       split_table = res$scores)
}

apply_elbow_post_corrections <- function(he, comp_df, K_init,
                                         h_elbow_redundant_jump_max = NULL,
                                         h_elbow_collapse_jump_ratio = Inf) {
  tab <- he$split_table
  if (is.null(tab) || nrow(tab) < 3L || is.na(he$K_opt)) return(he)

  # Overfit correction from offline_elbow_tuning.R:
  # K_init -> K_init - 1 is redundant when the first loss jump is tiny.
  if (!is.null(h_elbow_redundant_jump_max) &&
      is.finite(h_elbow_redundant_jump_max) &&
      he$K_opt == K_init - 1L) {
    idx_init <- which(tab$K == K_init)
    idx_target <- which(tab$K == K_init - 2L)

    if (length(idx_init) == 1L && length(idx_target) == 1L) {
      init_drop <- tab$nH2_jump[idx_init]
      if (is.finite(init_drop) &&
          init_drop <= h_elbow_redundant_jump_max) {
        he$K_opt <- tab$K[idx_target]
        he$best_idx <- comp_df$best_idx[idx_target]
        he$reason <- sprintf("overfit_one_step_corrected<=%.3g",
                             h_elbow_redundant_jump_max)
      }
    }
  }

  # Underfit correction from offline_elbow_tuning.R:
  # If K_init - 3 is selected only because the next drop collapses
  # catastrophically, keep the pre-collapse K_init - 2 solution.
  if (is.finite(h_elbow_collapse_jump_ratio) &&
      he$K_opt == K_init - 3L) {
    idx_prev <- which(tab$K == K_init - 2L)
    idx_sel <- which(tab$K == K_init - 3L)

    if (length(idx_prev) == 1L && length(idx_sel) == 1L) {
      prev_jump <- tab$nH2_jump[idx_prev]
      collapse_jump <- tab$nH2_jump[idx_sel]
      jump_ratio <- collapse_jump / pmax(prev_jump, 1e-12)

      if (is.finite(jump_ratio) &&
          jump_ratio >= h_elbow_collapse_jump_ratio) {
        he$K_opt <- tab$K[idx_prev]
        he$best_idx <- comp_df$best_idx[idx_prev]
        he$reason <- sprintf("underfit_collapse_corrected_ratio>=%.3g",
                             h_elbow_collapse_jump_ratio)
      }
    }
  }

  he
}

select_cached_pmhd_fit <- function(pmhd_cache, n_obs, K_init,
                                   h_elbow_cross_frac = 0.70,
                                   h_elbow_imp_min = 0.20,
                                   h_elbow_flat_ratio = 2.0,
                                   h_elbow_jump_frac = 0.35,
                                   h_elbow_jump_min = 1.0,
                                   h_elbow_redundant_jump_max = NULL,
                                   h_elbow_collapse_jump_ratio = Inf) {
  if (is.null(pmhd_cache$compressed_df) ||
      is.null(pmhd_cache$compressed_fits)) {
    return(list(fit = NULL, elbow = NULL))
  }

  he <- choose_H_elbow(pmhd_cache$compressed_df,
    cross_frac = h_elbow_cross_frac,
    imp_min = h_elbow_imp_min,
    flat_ratio = h_elbow_flat_ratio,
    n_obs = n_obs,
    jump_frac = h_elbow_jump_frac,
    jump_min = h_elbow_jump_min)
  he <- apply_elbow_post_corrections(
    he, pmhd_cache$compressed_df, K_init,
    h_elbow_redundant_jump_max = h_elbow_redundant_jump_max,
    h_elbow_collapse_jump_ratio = h_elbow_collapse_jump_ratio)

  fit_pos <- match(he$best_idx, pmhd_cache$compressed_df$best_idx)
  fit <- if (is.na(fit_pos)) NULL else pmhd_cache$compressed_fits[[fit_pos]]
  list(fit = fit, elbow = he)
}


# ============================================================
# Wrapper
# FIX: default elbow thresholds now match paper
# ============================================================

select_model_stability <- function(
    x, K_init = 5,
    lambda_grid = NULL,  # NULL = auto-generate coarse grid
    gamma = 3, nstart = 10,
    grid_n = 512, bw_adjust = 0.80,
    tol = 1e-5, max_iter = 100,
    delta_inner = 1e-8, delta_final = 1e-3,
    maxeval = 250,
    alpha_lower = 1.0, alpha_upper = 3.0,
    tau_lower = 0.01, tau_upper = 0.99,
    seed = 123,
    selection = c("H_elbow", "aic", "bic"),
    bic_params_per_K = 4L,
    bic_penalty_mult = 1.0,
    aic_penalty_mult = 1.0,
    K_min = 1L,
    h_elbow_cross_frac = 0.70,
    h_elbow_imp_min    = 0.20,
    h_elbow_flat_ratio = 2.0,
    h_elbow_jump_frac  = 0.35,
    h_elbow_jump_min   = 1.0,
    # Optional post-elbow corrections used by tuned Scenario 1.
    # Leave disabled by default so other scenarios keep the base selector.
    h_elbow_redundant_jump_max = NULL,
    h_elbow_collapse_jump_ratio = Inf,
    fixed_K_refit = TRUE,
    fixed_K_refit_nstart = min(nstart, 5),
    draw_paths = FALSE, plot_file = NULL,
    parallel = TRUE, ncores = 4) {

  selection <- match.arg(selection)
  if (!is.null(seed)) set.seed(seed)

  # --- Bandwidth: use bw.nrd0 (Silverman's rule) ---
  # bw.nrd0 internally uses min(sd, IQR/1.34) which is already 
  # moderately robust. The bandwidth itself is NOT the problem.
  bw0 <- bw.nrd0(x) * bw_adjust

  # --- Grid range: ROBUST to outliers ---
  # The original code used from=min(x)-3*bw0, to=max(x)+3*bw0.
  # Under Cauchy contamination, a single x=卤500 outlier stretches the 
  # grid to cover [-500, 500], spreading 512 grid points over 1000 units 
  # so the true mixture region [-8,8] gets only ~8 grid points.
  # Fix: anchor the grid to robust location/scale, ignoring outliers.
  robust_loc   <- median(x)
  robust_scale <- IQR(x) / 1.349  # consistent estimator of sigma
  if (!is.finite(robust_scale) || robust_scale < 1e-6)
    robust_scale <- mad(x, constant = 1.4826)
  grid_half <- max(4 * robust_scale, 3 * bw0, 6)  # at least 卤6
  grid_from <- robust_loc - grid_half - 3 * bw0
  grid_to   <- robust_loc + grid_half + 3 * bw0

  dens <- density(x, n = grid_n, bw = bw0,
                  from = grid_from, to = grid_to)
  grid <- dens$x; wq <- trapz_weights(grid)
  ghat <- normalize_grid_density(dens$y, wq)

  # --- Auto-generate coarse lambda grid if not provided ---
  # Stage 1 grid: SPARSE, covers a wide range [0, 3.0].
  # Goal: spend minimal points in the K=K_init plateau (where all
  # fits give the same overfitted answer), and let Stage 2 
  # adaptively refine around the K transitions.
  if (is.null(lambda_grid)) {
    lambda_grid <- unique(sort(c(0,
      exp(seq(log(1e-4), log(0.01), length.out = 3)),   # very sparse in K=K_init zone
      exp(seq(log(0.01), log(0.10), length.out = 5)),   # moderate near first transition
      exp(seq(log(0.10), log(1.00), length.out = 8)),   # moderate in main transition
      exp(seq(log(1.00), log(3.00), length.out = 4))))) # sparse at large 位
  }

  run_lambda <- function(i, lg = lambda_grid) {
    best_h2 <- Inf; best_f <- NULL
    for (s in seq_len(nstart)) {
      ai <- if (s == 1) 2.0 else
              runif(1, max(alpha_lower, 0.5), min(alpha_upper, 3.0))
      f <- tryCatch(
        pmhd_mcp_fit_one(
          x, K_init, lg[i], gamma, grid, wq, ghat,
          tol = tol, max_iter = max_iter,
          delta_inner = delta_inner, delta_final = delta_final,
          init = if (s == 1) NULL
                 else init_params_kmeans(x, K_init, alpha_init = ai),
          maxeval = maxeval,
          alpha_lower = alpha_lower, alpha_upper = alpha_upper,
          tau_lower = tau_lower, tau_upper = tau_upper),
        error = function(e) NULL)
      if (!is.null(f) && f$H2 < best_h2) {
        best_h2 <- f$H2; best_f <- f
      }
    }
    if (is.null(best_f))
      list(K = NA_integer_, H2 = Inf, H = Inf, fit = NULL)
    else
      list(K = best_f$Khat, H2 = best_f$H2,
           H = sqrt(best_f$H2), fit = best_f)
  }

  # --- STAGE 1: Run the provided lambda_grid ---
  if (parallel) {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(stop_parallel_backend(cl), add = TRUE)
    .efns <- c("pmhd_mcp_fit_one", "init_params_kmeans", "daepd",
               "trapz_weights", "normalize_grid_density", "H2_grid",
               "mcp_value", "mcp_prime")
    raw_results <- foreach(
      i = seq_along(lambda_grid),
      .packages = c("nloptr", "MASS"), .export = .efns
    ) %dopar% run_lambda(i)
  } else {
    raw_results <- lapply(seq_along(lambda_grid), run_lambda)
  }

  K_vec  <- sapply(raw_results, `[[`, "K")
  H2_vec <- sapply(raw_results, `[[`, "H2")
  H_vec  <- sapply(raw_results, `[[`, "H")

  if (all(is.na(K_vec))) {
    warning("select_model_stability: all fits failed"); return(NULL)
  }
  K_vec[is.na(K_vec)] <- max(K_vec, na.rm = TRUE)

  # --- STAGE 2: Adaptive refinement at K transition boundaries ---
  # Find lambda values where K changes, add denser grid points there.
  K_sorted <- K_vec[order(lambda_grid)]
  lam_sorted <- sort(lambda_grid)
  transitions <- which(diff(K_sorted) != 0)
  if (length(transitions) > 0) {
    refine_lams <- numeric(0)
    for (ti in transitions) {
      lam_lo <- lam_sorted[ti]
      lam_hi <- lam_sorted[min(ti + 1, length(lam_sorted))]
      if (lam_hi > lam_lo && lam_lo > 0) {
        new_lams <- exp(seq(log(max(lam_lo, 1e-6)),
                            log(lam_hi), length.out = 6))
        refine_lams <- c(refine_lams, new_lams[2:5])  # skip endpoints
      }
    }
    refine_lams <- setdiff(round(refine_lams, 8), round(lambda_grid, 8))
    if (length(refine_lams) > 0) {
      refine_results <- lapply(seq_along(refine_lams), function(i) {
        run_lambda(i, lg = refine_lams)
      })
      lambda_grid <- c(lambda_grid, refine_lams)
      raw_results <- c(raw_results, refine_results)
      K_vec  <- c(K_vec, sapply(refine_results, `[[`, "K"))
      H2_vec <- c(H2_vec, sapply(refine_results, `[[`, "H2"))
      H_vec  <- c(H_vec, sapply(refine_results, `[[`, "H"))
      # Re-sort by lambda
      ord <- order(lambda_grid)
      lambda_grid <- lambda_grid[ord]
      raw_results <- raw_results[ord]
      K_vec <- K_vec[ord]; H2_vec <- H2_vec[ord]; H_vec <- H_vec[ord]
    }
  }

  # --- STAGE 2b: Extend lambda upward if K at max(lambda) is still > 1 ---
  # If the largest lambda didn't push K down to at least 2 distinct 
  # levels, the grid doesn't cover the full transition. Extend upward.
  K_at_max <- K_vec[which.max(lambda_grid)]
  if (!is.na(K_at_max) && K_at_max > 1) {
    lam_max <- max(lambda_grid)
    extend_lams <- lam_max * c(1.5, 2.0, 3.0, 5.0)
    extend_results <- lapply(seq_along(extend_lams), function(i) {
      run_lambda(i, lg = extend_lams)
    })
    lambda_grid <- c(lambda_grid, extend_lams)
    raw_results <- c(raw_results, extend_results)
    K_vec  <- c(K_vec, sapply(extend_results, `[[`, "K"))
    H2_vec <- c(H2_vec, sapply(extend_results, `[[`, "H2"))
    H_vec  <- c(H_vec, sapply(extend_results, `[[`, "H"))
    ord <- order(lambda_grid)
    lambda_grid <- lambda_grid[ord]
    raw_results <- raw_results[ord]
    K_vec <- K_vec[ord]; H2_vec <- H2_vec[ord]; H_vec <- H_vec[ord]
  }

  K_vec[is.na(K_vec)] <- max(K_vec, na.rm = TRUE)
  K_mono <- K_vec
  for (i in 2:length(K_mono))
    K_mono[i] <- min(K_mono[i], K_mono[i - 1])

  n_obs <- length(x)
  df_vec <- sapply(raw_results, function(r) {
    if (is.null(r$fit) || length(r$fit$pi) == 0) return(Inf)
    r$fit$Khat * bic_params_per_K + (r$fit$Khat - 1)
  })
  nH2_vec   <- n_obs * H2_vec
  score_aic <- nH2_vec + aic_penalty_mult * 2 * df_vec
  score_bic <- nH2_vec + bic_penalty_mult * log(n_obs) * df_vec
  score_aic[!is.finite(score_aic)] <- Inf
  score_bic[!is.finite(score_bic)] <- Inf

  comp_df <- compress_lambda_path(lambda_grid, K_mono, H2_vec, H_vec,
                                  K_raw = K_vec)

  fixed_K_fit <- function(Km) {
    best_h2 <- Inf; best_fit <- NULL
    ntry <- max(1L, fixed_K_refit_nstart)
    for (s in seq_len(ntry)) {
      ai <- if (s == 1) 2.0 else
              runif(1, max(alpha_lower, 0.5), min(alpha_upper, 3.0))
      f <- tryCatch(
        pmhd_mcp_fit_one(
          x, Km, 0, gamma, grid, wq, ghat,
          tol = tol, max_iter = max_iter,
          delta_inner = delta_inner,
          delta_final = max(delta_final, 1e-2),
          init = if (s == 1) NULL
                 else init_params_kmeans(x, Km, alpha_init = ai),
          maxeval = maxeval,
          alpha_lower = alpha_lower, alpha_upper = alpha_upper,
          tau_lower = tau_lower, tau_upper = tau_upper),
        error = function(e) NULL)
      if (!is.null(f) && f$Khat == Km && f$H2 < best_h2) {
        best_h2 <- f$H2; best_fit <- f
      }
    }
    if (is.null(best_fit)) return(NULL)
    list(K = Km, H2 = best_h2, H = sqrt(best_h2), fit = best_fit)
  }

  insert_fixed_K <- function(Km) {
    if (!isTRUE(fixed_K_refit) || Km %in% comp_df$K) return(FALSE)
    gap <- fixed_K_fit(Km)
    if (is.null(gap)) return(FALSE)
    new_idx <- length(raw_results) + 1L
    raw_results[[new_idx]] <<- gap
    comp_df <<- rbind(comp_df, data.frame(
      K = Km, W = 1L, lam_range = 0,
      best_idx = new_idx,
      lambda = 0, H2_star = gap$H2,
      H_star = gap$H,
      logH_star = log(max(gap$H, 1e-12))))
    TRUE
  }

  if (isTRUE(fixed_K_refit) && nrow(comp_df) == 1L &&
      comp_df$K[1] == K_init && K_init > K_min) {
    for (Km in seq(K_init - 1L, K_min, by = -1L)) {
      insert_fixed_K(Km)
    }
    comp_df <- comp_df[order(comp_df$K, decreasing = TRUE), , drop = FALSE]
    rownames(comp_df) <- NULL
  }

  # --- Targeted gap-fill for K values skipped by biggest 螖K jump ---
  # If the compressed path has a jump 螖K 鈮?2 (e.g., K=4鈫扠=2), the 
  # intermediate K values (K=3) never appeared in the MCP path.
  # Fit these missing K values explicitly with lambda=0 and insert 
  # them into comp_df so the selection rule can consider them.
  if (nrow(comp_df) >= 2) {
    Kvals_cd <- comp_df$K  # decreasing
    deltas_cd <- c(0, Kvals_cd[-length(Kvals_cd)] - Kvals_cd[-1])
    big_jumps <- which(deltas_cd >= 2)
    
    for (ji in big_jumps) {
      K_above <- Kvals_cd[ji - 1]  # K before the jump
      K_below <- Kvals_cd[ji]       # K after the jump
      missing_in_jump <- (K_below + 1):(K_above - 1)  # K values skipped
      
      for (Km in missing_in_jump) {
        insert_fixed_K(Km)
      }
    }
    comp_df <- comp_df[order(comp_df$K, decreasing = TRUE), , drop = FALSE]
    rownames(comp_df) <- NULL
  }

  h_elb   <- choose_H_elbow(comp_df,
                             cross_frac = h_elbow_cross_frac,
                             imp_min    = h_elbow_imp_min,
                             flat_ratio = h_elbow_flat_ratio,
                             n_obs = n_obs,
                             jump_frac = h_elbow_jump_frac,
                             jump_min  = h_elbow_jump_min)
  h_elb <- apply_elbow_post_corrections(
    h_elb, comp_df, K_init,
    h_elbow_redundant_jump_max = h_elbow_redundant_jump_max,
    h_elbow_collapse_jump_ratio = h_elbow_collapse_jump_ratio)
  idx_He <- h_elb$best_idx

  pick <- function(score) {
    v <- K_mono >= K_min
    if (any(v) && any(is.finite(score[v]))) {
      vi <- which(v); vi[which.min(score[vi])]
    } else which.min(score)
  }
  idx_aic <- pick(score_aic)
  idx_bic <- pick(score_bic)

  best_idx <- switch(selection,
                     H_elbow = idx_He, aic = idx_aic, bic = idx_bic)

  path_df <- data.frame(
    lambda = lambda_grid,
    log10_lambda = log10(pmax(lambda_grid, 1e-12)),
    K_raw = K_vec, K_mono = K_mono,
    H2 = H2_vec, H = H_vec, nH2 = nH2_vec,
    df = df_vec, score_aic = score_aic, score_bic = score_bic)

  # Safe accessor: best_idx may point to a gap-filled entry beyond
  # the original lambda_grid / K_mono length.
  safe_K <- function(idx) {
    if (is.na(idx) || idx > length(K_mono)) {
      # Gap-filled entry: get K from raw_results or comp_df
      r <- raw_results[[idx]]
      if (!is.null(r)) r$K else NA_integer_
    } else K_mono[idx]
  }
  safe_lambda <- function(idx) {
    if (is.na(idx) || idx > length(lambda_grid)) NA_real_
    else lambda_grid[idx]
  }
  safe_fit <- function(idx) {
    if (is.na(idx) || idx > length(raw_results)) NULL
    else raw_results[[idx]]$fit
  }
  compressed_fits <- lapply(comp_df$best_idx, safe_fit)
  names(compressed_fits) <- as.character(comp_df$K)

  list(fit_opt = safe_fit(best_idx),
       K_opt = h_elb$K_opt,  # use elbow's own K, not K_mono
       lambda_opt = safe_lambda(best_idx),
       best_idx = best_idx,
       K_H_elbow = h_elb$K_opt,
       lambda_H_elbow = safe_lambda(idx_He),
       K_aic = safe_K(idx_aic), lambda_aic = safe_lambda(idx_aic),
       K_bic = safe_K(idx_bic), lambda_bic = safe_lambda(idx_bic),
       path_df = path_df, compressed_df = comp_df,
       compressed_fits = compressed_fits,
       h_elbow_detail = h_elb)
}


# ============================================================
# Helpers
# ============================================================

# align_labels <- function(est_params, true_params) {
#   K <- nrow(true_params)
#   dist_mat <- matrix(0, K, K)
#   for (i in 1:K) for (j in 1:K)
#     dist_mat[i, j] <- sum((est_params[i, ] - true_params[j, ])^2)
#   clue::solve_LSAP(dist_mat)
# }

align_labels <- function(est_params, true_params) {
  K <- nrow(true_params)
  dist_mat <- matrix(0, K, K)
  for (i in 1:K) for (j in 1:K)
    dist_mat[i, j] <- sum((est_params[i, ] - true_params[j, ])^2)
  order(as.integer(clue::solve_LSAP(dist_mat)))
}

compute_mse_aligned <- function(est, true_pi, true_mu, true_sigma,
                                true_alpha, true_tau) {
  K0 <- length(true_mu)
  if (est$Khat != K0) return(list(mse = NA, bias = NA,
    mse_pi = NA, mse_mu = NA, mse_sigma = NA,
    mse_alpha = NA, mse_tau = NA))
  true_mat <- cbind(true_pi, true_mu, true_sigma, true_alpha, true_tau)
  est_mat  <- cbind(est$pi, est$mu, est$sigma, est$alpha, est$tau)
  perm <- align_labels(est_mat, true_mat)
  est_aligned <- est_mat[perm, ]
  diffs <- est_aligned - true_mat
  list(mse  = mean(diffs^2),
       bias = colMeans(diffs),
       mse_pi    = mean(diffs[, 1]^2),
       mse_mu    = mean(diffs[, 2]^2),
       mse_sigma = mean(diffs[, 3]^2),
       mse_alpha = mean(diffs[, 4]^2),
       mse_tau   = mean(diffs[, 5]^2),
       perm = perm)
}

aepd_component_moments <- function(mu, sigma, alpha, tau) {
  alpha <- pmax(alpha, 1e-10)
  tau <- pmin(pmax(tau, 1e-8), 1 - 1e-8)
  ev1 <- gamma(2 / alpha) / gamma(1 / alpha)
  ev2 <- gamma(3 / alpha) / gamma(1 / alpha)
  shift <- sigma * ev1 * ((1 - tau) / tau - tau / (1 - tau))
  raw2 <- sigma^2 * ev2 * ((1 - tau) / tau^2 + tau / (1 - tau)^2)
  var <- pmax(raw2 - shift^2, 0)
  list(mean = mu + shift, sd = sqrt(var))
}

# MSE for GM-MLE against the first two moments of the true AEPD components.
compute_mse_gm <- function(mc, true_pi, true_mu, true_sigma_aepd,
                           true_alpha = rep(2, length(true_mu)),
                           true_tau = rep(0.5, length(true_mu))) {
  K0 <- length(true_mu)
  if (is.null(mc) || mc$G != K0) return(list(mse = NA, mse_pi = NA,
    mse_mu = NA, mse_sd = NA))
  pro <- mc$parameters$pro
  if (is.null(pro)) pro <- 1
  sr <- mc$parameters$variance$sigmasq
  if (is.null(sr)) sr <- mc$parameters$variance$sigma2
  sd_mc <- sqrt(if (length(sr) == 1L) rep(sr, mc$G) else sr)

  mom <- aepd_component_moments(true_mu, true_sigma_aepd,
                                true_alpha, true_tau)
  true_mat <- cbind(true_pi, mom$mean, mom$sd)
  est_mat  <- cbind(pro, mc$parameters$mean, sd_mc)
  perm <- align_labels(est_mat, true_mat)
  est_aligned <- est_mat[perm, ]
  diffs <- est_aligned - true_mat
  list(mse    = mean(diffs^2),
       mse_pi = mean(diffs[, 1]^2),
       mse_mu = mean(diffs[, 2]^2),
       mse_sd = mean(diffs[, 3]^2))
}

true_density_grid <- function(grid, pi, mu, sigma, alpha, tau) {
  wq <- trapz_weights(grid)
  f  <- rowSums(sapply(seq_along(pi), function(k)
    pi[k] * daepd(grid, mu[k], sigma[k], alpha[k], tau[k])))
  f <- pmax(f, 1e-300); f / sum(f * wq)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

# --- AEPD density on grid ---
aepd_density_grid <- function(grid, fit) {
  if (is.null(fit)) return(rep(NA, length(grid)))
  f <- rowSums(sapply(seq_along(fit$pi), function(k)
    fit$pi[k] * daepd(grid, fit$mu[k], fit$sigma[k],
                      fit$alpha[k], fit$tau[k])))
  pmax(f, 1e-300)
}

# --- GM density on grid ---
gm_density_grid <- function(grid, mc) {
  if (is.null(mc)) return(rep(NA, length(grid)))
  pro <- mc$parameters$pro
  if (is.null(pro)) pro <- 1
  sr <- mc$parameters$variance$sigmasq
  if (is.null(sr)) sr <- mc$parameters$variance$sigma2
  sig <- sqrt(if (length(sr) == 1L) rep(sr, mc$G) else sr)
  f <- rowSums(vapply(seq_len(mc$G), function(k)
    pro[k] * dnorm(grid, mc$parameters$mean[k], sig[k]),
    numeric(length(grid))))
  pmax(f, 1e-300)
}


# ============================================================
# Random generation
# ============================================================

raepd <- function(n, mu, sigma, alpha, tau) {
  u <- runif(n); s <- ifelse(u < (1 - tau), 1L, -1L)
  g <- rgamma(n, shape = 1 / alpha, rate = 1)
  v <- g^(1 / alpha)
  ts <- ifelse(s == 1L, tau, 1 - tau)
  mu + s * sigma * v / ts
}

rmix_aepd <- function(n, pi, mu, sigma, alpha, tau) {
  K <- length(pi); ni <- as.vector(rmultinom(1, n, pi))
  x <- numeric(n); idx <- 0
  for (k in seq_len(K)) {
    if (ni[k] == 0) next
    x[(idx + 1):(idx + ni[k])] <-
      raepd(ni[k], mu[k], sigma[k], alpha[k], tau[k])
    idx <- idx + ni[k]
  }
  sample(x)
}

rskt <- function(n, xi = 0, omega = 1, alpha_skew = 0, nu = 3) {
  delta <- alpha_skew / sqrt(1 + alpha_skew^2)
  u  <- rgamma(n, shape = nu / 2, rate = nu / 2)
  z0 <- rnorm(n); z1 <- rnorm(n)
  y  <- delta * abs(z0) / sqrt(u) +
        sqrt(1 - delta^2) * z1 / sqrt(u)
  xi + omega * y
}

dskt <- function(x, xi = 0, omega = 1, alpha_skew = 0, nu = 3) {
  z <- (x - xi) / omega
  d <- 2 / omega * dt(z, df = nu) *
    pt(alpha_skew * z * sqrt((nu + 1) / (nu + z^2)),
       df = nu + 1)
  pmax(d, 1e-300)
}


# ============================================================

# Source one or more scenario files after this file, for example:
#   source("scenario5.R")
#   res5 <- run_scenario5(B=100, parallel=TRUE, ncores=6)
