# ============================================================
# galaxies_final_analysis.R
# ============================================================
# Self-contained final empirical analysis for the galaxies data.
# Runs PMHD-MCP order selection, fixed-order PMHD-MCP-R refit,
# Gaussian-mixture comparison, and writes all final tables/figures.
# Default output directory: ./K6 next to this script.
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b
script_path <- {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    sub("^--file=", "", file_arg[1L])
  } else {
    sys.frame(1)$ofile %||% getwd()
  }
}
galaxies_dir <- dirname(normalizePath(script_path, winslash = "/",
                                      mustWork = FALSE))
galaxies_dir <- Sys.getenv("GALAXIES_DIR", unset = galaxies_dir)
galaxies_dir <- normalizePath(galaxies_dir, winslash = "/", mustWork = FALSE)

# ---- Core PMHD-MCP functions --------------------------------
# ============================================================
# Galaxies Dataset Analysis with PMHD-MCP AEPD
# ------------------------------------------------------------
# Synced to the H_elbow / H-AIC / H-BIC framework used in the
# updated simulation code.
#
# Main changes:
#   1. Scan the lambda path ONCE and cache the best fit at each lambda
#   2. Compute three Hellinger-based selectors simultaneously:
#        - H_elbow  (default)
#        - H-AIC    = n * H^2 + aic_penalty_mult * 2 * df
#        - H-BIC    = n * H^2 + bic_penalty_mult * log(n) * df
#   3. Draw selector path plots
#
# Dataset : MASS::galaxies (n = 82 recession velocities in km/s)
# Competitor: Gaussian mixture via mclust (BIC)
# ============================================================

suppressPackageStartupMessages({
  library(MASS)
  library(nloptr)
  library(mclust)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------------------------------------
# Core utilities
# ------------------------------------------------------------
trapz_weights <- function(x) {
  n <- length(x); dx <- diff(x); w <- numeric(n)
  w[1] <- dx[1] / 2; w[n] <- dx[n - 1] / 2
  if (n > 2) w[2:(n - 1)] <- (dx[1:(n - 2)] + dx[2:(n - 1)]) / 2
  w
}

normalize_grid_density <- function(y, w, eps = 1e-15) {
  y <- pmax(y, 0); z <- sum(y * w)
  if (!is.finite(z) || z <= eps) return(rep(1 / length(y), length(y)))
  y / z
}

H2_grid <- function(p, q, w) {
  p <- pmax(p, 0); q <- pmax(q, 0)
  val <- 1 - sum(sqrt(p * q) * w)
  2 * pmax(val, 0)
}

mcp_value <- function(pi, lambda, gamma)
  sum(ifelse(pi <= gamma * lambda,
             lambda * pi - pi^2 / (2 * gamma),
             gamma * lambda^2 / 2))

mcp_prime <- function(pi, lambda, gamma) pmax(lambda - pi / gamma, 0)

# ------------------------------------------------------------
# AEPD density and initialization
# ------------------------------------------------------------
daepd <- function(x, mu, sigma, alpha, tau) {
  sigma <- pmax(sigma, 1e-10); alpha <- pmax(alpha, 1e-10)
  tau   <- pmin(pmax(tau, 1e-8), 1 - 1e-8)
  C <- alpha * tau * (1 - tau) / (base::gamma(1 / alpha) * sigma)
  z <- abs(x - mu)^alpha / sigma^alpha
  d <- C * exp(-z * ifelse(x < mu, (1 - tau)^alpha, tau^alpha))
  pmax(d, 1e-300)
}

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

# ------------------------------------------------------------
# PMHD-MCP AEPD core fitter
# ------------------------------------------------------------
pmhd_mcp_fit_one <- function(x, K_init, lambda, gamma = 3, grid, wq, ghat,
                             tol = 1e-6, max_iter = 200,
                             delta_inner = 1e-8, delta_final = 1e-3,
                             init = NULL, maxeval = 400,
                             alpha_lower = 1.0, alpha_upper = 8.5,
                             tau_lower = 0.01, tau_upper = 0.99) {

  sigma_lower <- 1e-3 * sd(x); sigma_upper <- 3 * diff(range(x))
  log_sl <- log(sigma_lower); log_su <- log(sigma_upper)
  log_al <- log(alpha_lower); log_au <- log(alpha_upper)

  if (is.null(init)) init <- init_params_kmeans(x, K_init)
  pi_k <- pmax(init$pi, delta_inner); pi_k <- pi_k / sum(pi_k)
  mu_k <- init$mu; sigma_k <- init$sigma
  alpha_k <- init$alpha; tau_k <- init$tau
  obj_old <- Inf; H2_now <- NA_real_

  for (iter in seq_len(max_iter)) {
    # E-step
    f_old <- matrix(0, length(grid), K_init)
    for (k in seq_len(K_init))
      f_old[, k] <- daepd(grid, mu_k[k], sigma_k[k], alpha_k[k], tau_k[k])
    f_mix_old <- pmax(drop(f_old %*% pi_k), 1e-300)
    R_old <- sweep(sweep(f_old, 2, pi_k, "*"), 1, f_mix_old, "/")
    c_k <- colSums(ghat * R_old * wq)

    # M-step for theta
    for (k in which(pi_k >= delta_inner)) {
      gtilde_k <- (ghat * R_old[, k]) / pmax(c_k[k], 1e-12)

      obj_theta <- function(par) {
        mu <- par[1]; sig <- exp(par[2])
        alp <- exp(par[3]); t <- plogis(par[4])
        1 - sum(sqrt(gtilde_k * daepd(grid, mu, sig, alp, t)) * wq)
      }

      res <- tryCatch(
        nloptr(c(mu_k[k], log(sigma_k[k]), log(alpha_k[k]), qlogis(tau_k[k])),
               obj_theta,
               lb = c(-Inf, log_sl, log_al, qlogis(tau_lower)),
               ub = c(Inf,  log_su, log_au, qlogis(tau_upper)),
               opts = list(algorithm = "NLOPT_LN_BOBYQA",
                           xtol_rel = 1e-8, maxeval = maxeval)),
        error = function(e) NULL)
      if (!is.null(res)) {
        mu_k[k]    <- res$solution[1]
        sigma_k[k] <- exp(res$solution[2])
        alpha_k[k] <- exp(res$solution[3])
        tau_k[k]   <- plogis(res$solution[4])
      }
    }

    # M-step for pi
    f_new <- matrix(0, length(grid), K_init)
    for (k in seq_len(K_init))
      f_new[, k] <- daepd(grid, mu_k[k], sigma_k[k], alpha_k[k], tau_k[k])
    B_k <- 2 * colSums(sqrt(pmax(ghat * R_old * f_new, 0)) * wq)
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
    if (sum(pi_k) <= 0) pi_k <- rep(1 / K_init, K_init) else pi_k <- pi_k / sum(pi_k)

    # Convergence check
    H2_now <- H2_grid(ghat, pmax(drop(f_new %*% pi_k), 1e-300), wq)
    obj_now <- H2_now + mcp_value(pi_k, lambda, gamma)
    if (abs(obj_old - obj_now) < tol) break
    obj_old <- obj_now
  }

  keep <- which(pi_k >= delta_final)
  if (length(keep) == 0L) keep <- which.max(pi_k)
  pi_keep <- pi_k[keep] / sum(pi_k[keep])
  f_keep <- matrix(0, length(grid), length(keep))
  for (j in seq_along(keep)) {
    k <- keep[j]
    f_keep[, j] <- daepd(grid, mu_k[k], sigma_k[k], alpha_k[k], tau_k[k])
  }
  H2_final <- H2_grid(ghat, pmax(drop(f_keep %*% pi_keep), 1e-300), wq)

  list(pi = pi_keep,
       mu = mu_k[keep], sigma = sigma_k[keep],
       alpha = alpha_k[keep], tau = tau_k[keep],
       Khat = length(keep), H2 = H2_now, H2_final = H2_final,
       n_iter = iter)
}

# ------------------------------------------------------------
# H_elbow path processing (same framework as simulation)
# ------------------------------------------------------------
weighted_lm_sse <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & (w > 0)
  x <- x[ok]; y <- y[ok]; w <- w[ok]
  if (length(y) <= 1L) return(0)
  fit <- tryCatch(stats::lm.wfit(cbind(1, x), y, w = w), error = function(e) NULL)
  if (is.null(fit)) return(Inf)
  sum(w * fit$residuals^2)
}

compress_lambda_path <- function(lambda_grid, K_mono, H2_vec, H_vec) {
  K_levels <- sort(unique(K_mono[is.finite(K_mono)]), decreasing = TRUE)
  if (length(K_levels) == 0L) {
    return(data.frame(K = numeric(0), W = integer(0),
                      lam_range = numeric(0),
                      best_idx = integer(0),
                      lambda = numeric(0), H2_star = numeric(0),
                      H_star = numeric(0), logH_star = numeric(0)))
  }
  lambda_pos <- lambda_grid[is.finite(lambda_grid) & lambda_grid > 0]
  lambda_floor <- if (length(lambda_pos) > 0L) min(lambda_pos) / 2 else 1e-12
  log_lambda <- log(pmax(lambda_grid, lambda_floor))
  rows <- lapply(K_levels, function(Kv) {
    idx <- which(K_mono == Kv & is.finite(H2_vec))
    best_idx <- idx[which.min(H2_vec[idx])]
    log_lams_at_K <- log_lambda[idx]
    data.frame(
      K = Kv,
      W = length(idx),
      lam_range = max(log_lams_at_K) - min(log_lams_at_K),
      best_idx = best_idx,
      lambda = lambda_grid[best_idx],
      H2_star = H2_vec[best_idx],
      H_star = H_vec[best_idx],
      logH_star = log(pmax(H_vec[best_idx], 1e-12))
    )
  })
  do.call(rbind, rows)
}

choose_K_from_path <- function(K_mono, compressed_df, n_obs = NULL,
                               jump_frac = 0.03, jump_min = 0) {
  if (is.null(compressed_df) || nrow(compressed_df) == 0L)
    return(list(K_opt = NA_integer_, best_idx = NA_integer_,
                reason = "empty_path", scores = NULL))

  compressed_df <- compressed_df[is.finite(compressed_df$H_star), ,
                                 drop = FALSE]
  J <- nrow(compressed_df)
  if (J == 0L)
    return(list(K_opt = NA_integer_, best_idx = NA_integer_,
                reason = "empty_path", scores = NULL))
  if (J == 1L)
    return(list(K_opt = compressed_df$K[1L],
                best_idx = compressed_df$best_idx[1L],
                reason = "single_K", scores = NULL))

  Kvals <- compressed_df$K
  Wvals <- compressed_df$W
  Rvals <- compressed_df$lam_range
  H2vals <- compressed_df$H2_star
  nn <- if (!is.null(n_obs)) n_obs else 1
  nH2vals <- nn * H2vals

  delta_K <- c(0, Kvals[-J] - Kvals[-1L])
  nH2_last <- if (Kvals[J] > 1L) nH2vals[J] else 0
  nH2_jump <- c(pmax(nH2vals[-1L] - nH2vals[-J], 0), nH2_last)

  R_max <- max(Rvals, na.rm = TRUE)
  R_norm <- if (R_max > 1e-10) Rvals / R_max else rep(0, J)
  nH2_max <- max(nH2_jump, na.rm = TRUE)
  nH2_norm <- if (nH2_max > 1e-10) nH2_jump / nH2_max else rep(0, J)
  scores <- 3 * nH2_norm + R_norm
  scores[1L] <- -Inf
  if (Kvals[J] <= 1L) scores[J] <- -Inf

  jump_cut <- max(jump_min, jump_frac * nH2_max)
  candidates <- which(nH2_jump >= jump_cut)
  candidates <- candidates[candidates > 1L & Kvals[candidates] > 1L]
  if (length(candidates) > 0L) {
    pick <- candidates[1L]
    reason <- "first_substantial_nH2_jump"
  } else {
    finite_scores <- is.finite(scores)
    if (any(finite_scores)) {
      best_score <- max(scores[finite_scores], na.rm = TRUE)
      tied <- which(abs(scores - best_score) < 1e-10)
      pick <- tied[1L]
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
                           imp_min = 0.20,
                           flat_ratio = 2.0,
                           n_obs = NULL,
                           jump_frac = 0.03,
                           jump_min = 0) {
  res <- choose_K_from_path(NULL, compressed_df, n_obs = n_obs,
                            jump_frac = jump_frac, jump_min = jump_min)
  list(best_idx = res$best_idx,
       K_opt = res$K_opt,
       split_idx = NA_integer_,
       fallback = (res$reason == "single_K" || res$reason == "empty_path"),
       reason = res$reason,
       split_table = res$scores)
}

plot_selector_paths <- function(sel, file = NULL,
                                include_aic_bic = FALSE,
                                width = if (include_aic_bic) 10 else 10,
                                height = if (include_aic_bic) 8 else 4.8) {
  stopifnot(!is.null(sel$path_df), !is.null(sel$compressed_df))

  if (!is.null(file)) grDevices::pdf(file, width = width, height = height)
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    if (!is.null(file)) grDevices::dev.off()
  }, add = TRUE)

  path <- sel$path_df
  comp <- sel$compressed_df
  lx <- log10(pmax(path$lambda, 1e-12))

  if (isTRUE(include_aic_bic)) {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 4))
  } else {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 4))
  }

  # Panel 1: Hellinger path + K path + selected lambdas
  plot(lx, path$nH2, type = "b", pch = 16, col = "black",
       xlab = expression(log[10](lambda)),
       ylab = expression(n %.% H^2),
       main = "Hellinger loss path")
  abline(v = log10(sel$lambda_H_elbow), col = "firebrick", lty = 2, lwd = 2)
  if (isTRUE(include_aic_bic)) {
    abline(v = log10(sel$lambda_aic), col = "steelblue", lty = 3, lwd = 2)
    abline(v = log10(sel$lambda_bic), col = "forestgreen", lty = 4, lwd = 2)
  }
  par(new = TRUE)
  plot(lx, path$K_mono, type = "s", lwd = 1.5, col = "gray50",
       axes = FALSE, xlab = "", ylab = "")
  axis(4, col.axis = "gray50", col = "gray50")
  mtext("Monotone K", side = 4, line = 2.5, col = "gray50")

  if (isTRUE(include_aic_bic)) {
    legend("bottomleft", inset = c(0, 0.1),
           legend = c("n H^2", "K path", "H_elbow", "H-AIC", "H-BIC"),
           col = c("black", "gray50", "firebrick", "steelblue", "forestgreen"),
           lty = c(1, 1, 2, 3, 4),
           pch = c(16, NA, NA, NA, NA),
           bty = "n", cex = 0.8)
  } else {
    legend("bottomleft", inset = c(0, 0.1),
           legend = c("n H^2", "K path", "H_elbow"),
           col = c("black", "gray50", "firebrick"),
           lty = c(1, 1, 2),
           pch = c(16, NA, NA),
           bty = "n", cex = 0.8)
  }

  # Panel 2: H_elbow compressed path
  plot(comp$K, comp$logH_star, type = "b", pch = 16, lwd = 2,
       cex = 0.8 + 1.2 * comp$W / max(comp$W),
       xlab = "Compressed order K",
       ylab = expression(log(H[K]^"*")),
       main = "H_elbow compressed path")
  abline(v = sel$K_H_elbow, col = "firebrick", lty = 2, lwd = 2)

  if (isTRUE(include_aic_bic)) {
    # Panel 3: H-AIC
    plot(lx, path$score_aic, type = "b", pch = 16, col = "steelblue",
         xlab = expression(log[10](lambda)), ylab = "Score",
         main = "H-AIC path")
    abline(v = log10(sel$lambda_aic), col = "steelblue", lty = 2, lwd = 2)

    # Panel 4: H-BIC
    plot(lx, path$score_bic, type = "b", pch = 16, col = "forestgreen",
         xlab = expression(log[10](lambda)), ylab = "Score",
         main = "H-BIC path")
    abline(v = log10(sel$lambda_bic), col = "forestgreen", lty = 2, lwd = 2)
  }

  invisible(sel)
}


# ------------------------------------------------------------
# Single path scan + all selectors from Hellinger loss
# ------------------------------------------------------------
fit_pmhd_aepd_galaxies <- function(x,
                                   K_max = 8,
                                   lambda_grid = c(0, exp(seq(log(1e-4), log(0.5), length.out = 60))),
                                   gamma = 3,
                                   nstart = 20,
                                   grid_n = 1024,
                                   bw_adjust = 0.80,
                                   tol = 1e-6,
                                   max_iter = 200,
                                   delta_inner = 1e-8,
                                   delta_final = NULL,
                                   maxeval = 400,
                                   alpha_lower = 1.0,
                                   alpha_upper = 8.5,
                                   selector = c("H_elbow", "aic", "bic"),
                                   aic_penalty_mult = 1.0,
                                   bic_penalty_mult = 1.0,
                                   bic_params_per_K = 4L,
                                   h_elbow_jump_frac = 0.03,
                                   h_elbow_jump_min = 0.0,
                                   seed = 20250420,
                                   draw_paths = FALSE,
                                   plot_file = NULL,
                                   verbose = TRUE) {

  selector <- match.arg(selector)
  if (is.null(delta_final)) delta_final <- max(1e-3, 2 / length(x))
  if (!is.null(seed)) set.seed(seed)

  bw0  <- bw.nrd0(x) * bw_adjust
  dens <- density(x, n = grid_n, bw = bw0,
                  from = min(x) - 5 * bw0,
                  to   = max(x) + 5 * bw0)
  grid <- dens$x
  wq   <- trapz_weights(grid)
  ghat <- normalize_grid_density(dens$y, wq)

  if (verbose) {
    cat(sprintf("Using delta_final = %.4f (min effective n per component = %.2f)\n",
                delta_final, delta_final * length(x)))
    cat(sprintf("Fitting PMHD-MCP AEPD on n=%d, grid [%g, %g], bw=%.4f\n",
                length(x), min(grid), max(grid), bw0))
  }

  fits_by_lambda <- vector("list", length(lambda_grid))
  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    best_fit <- NULL; best_h2 <- Inf
    for (s in seq_len(nstart)) {
      init_i <- if (s == 1) NULL else init_params_kmeans(x, K_max)
      f <- tryCatch(
        pmhd_mcp_fit_one(x, K_max, lam, gamma, grid, wq, ghat,
                         tol = tol, max_iter = max_iter,
                         delta_inner = delta_inner,
                         delta_final = delta_final,
                         init = init_i, maxeval = maxeval,
                         alpha_lower = alpha_lower,
                         alpha_upper = alpha_upper),
        error = function(e) NULL)
      if (!is.null(f) && is.finite(f$H2) && f$H2 < best_h2) {
        best_h2 <- f$H2; best_fit <- f
      }
    }
    fits_by_lambda[[i]] <- best_fit
  }

  K_vec  <- sapply(fits_by_lambda, function(f) if (is.null(f)) NA_integer_ else f$Khat)
  H2_vec <- sapply(fits_by_lambda, function(f) if (is.null(f)) Inf else f$H2)
  H_vec  <- sqrt(H2_vec)

  if (all(is.na(K_vec))) stop("All lambda fits failed.")

  K_vec[is.na(K_vec)] <- max(K_vec, na.rm = TRUE)
  K_mono <- K_vec
  for (i in 2:length(K_mono)) K_mono[i] <- min(K_mono[i], K_mono[i - 1])

  n_obs <- length(x)
  df_vec <- sapply(fits_by_lambda, function(fit) {
    if (is.null(fit) || length(fit$pi) == 0) return(Inf)
    K <- fit$Khat
    K * bic_params_per_K + (K - 1)
  })

  nH2_vec <- n_obs * H2_vec
  score_aic <- nH2_vec + aic_penalty_mult * 2 * df_vec
  score_bic <- nH2_vec + bic_penalty_mult * log(n_obs) * df_vec
  score_aic[!is.finite(score_aic)] <- Inf
  score_bic[!is.finite(score_bic)] <- Inf

  compressed_df <- compress_lambda_path(lambda_grid, K_mono, H2_vec, H_vec)
  h_elbow <- choose_H_elbow(
    compressed_df,
    n_obs = n_obs,
    jump_frac = h_elbow_jump_frac,
    jump_min = h_elbow_jump_min
  )
  idx_H_elbow <- h_elbow$best_idx
  idx_aic <- which.min(score_aic)
  idx_bic <- which.min(score_bic)

  best_idx <- if (selector == "H_elbow") idx_H_elbow else if (selector == "aic") idx_aic else idx_bic

  path_df <- data.frame(
    idx = seq_along(lambda_grid),
    lambda = lambda_grid,
    log10_lambda = log10(pmax(lambda_grid, 1e-12)),
    K_raw = K_vec,
    K_mono = K_mono,
    H2 = H2_vec,
    H = H_vec,
    nH2 = nH2_vec,
    df = df_vec,
    score_aic = score_aic,
    score_bic = score_bic
  )

  best_fit_obj <- fits_by_lambda[[best_idx]]
  fit_H_elbow_obj <- fits_by_lambda[[idx_H_elbow]]
  fit_aic_obj <- fits_by_lambda[[idx_aic]]
  fit_bic_obj <- fits_by_lambda[[idx_bic]]

  out <- list(
    best_fit = fits_by_lambda[[best_idx]],
    best_lambda = lambda_grid[best_idx],
    best_idx = best_idx,
    best_K = best_fit_obj$Khat,
    selector = selector,

    fit_H_elbow = fits_by_lambda[[idx_H_elbow]],
    K_H_elbow = fit_H_elbow_obj$Khat,
    lambda_H_elbow = lambda_grid[idx_H_elbow],
    idx_H_elbow = idx_H_elbow,
    H_elbow_details = h_elbow,

    fit_aic = fits_by_lambda[[idx_aic]],
    K_aic = fit_aic_obj$Khat,
    lambda_aic = lambda_grid[idx_aic],
    idx_aic = idx_aic,

    fit_bic = fits_by_lambda[[idx_bic]],
    K_bic = fit_bic_obj$Khat,
    lambda_bic = lambda_grid[idx_bic],
    idx_bic = idx_bic,

    path_df = path_df,
    compressed_df = compressed_df,
    fits_by_lambda = fits_by_lambda,
    lambda_grid = lambda_grid,
    kde = list(grid = grid, wq = wq, ghat = ghat, bw = bw0),
    config = list(aic_penalty_mult = aic_penalty_mult,
                  bic_penalty_mult = bic_penalty_mult,
                  h_elbow_jump_frac = h_elbow_jump_frac,
                  h_elbow_jump_min = h_elbow_jump_min)
  )

  if (draw_paths) plot_selector_paths(out, file = plot_file)
  out
}

# ------------------------------------------------------------
# Gaussian mixture competitor
# ------------------------------------------------------------
fit_gm_mclust <- function(x, Kmax = 8) {
  mc <- Mclust(x, G = 1:Kmax, verbose = FALSE)
  s2 <- mc$parameters$variance$sigmasq
  if (length(s2) == 1L) s2 <- rep(s2, mc$G)
  sig <- sqrt(s2)
  list(Khat = mc$G, pi = mc$parameters$pro, mu = mc$parameters$mean,
       sigma = sig, model = mc)
}

eval_aepd_mixture <- function(grid, pi, mu, sigma, alpha, tau) {
  rowSums(sapply(seq_along(pi), function(k)
    pi[k] * daepd(grid, mu[k], sigma[k], alpha[k], tau[k])))
}

eval_gm <- function(grid, pi, mu, sigma) {
  rowSums(sapply(seq_along(pi), function(k)
    pi[k] * dnorm(grid, mu[k], sigma[k])))
}

describe_shape <- function(alpha, tau) {
  tail <- if (alpha > 2.1) "platykurtic" else if (alpha < 1.9) "leptokurtic" else "near-Gaussian"
  skew <- if (tau > 0.53) "right-skewed" else if (tau < 0.47) "left-skewed" else "near-symmetric"
  paste0(tail, ", ", skew)
}

# ------------------------------------------------------------
# Optional stability diagnostic under the chosen selector
# ------------------------------------------------------------
stability_check <- function(x, K_max = 8,
                            base_lambda_grid = c(0, exp(seq(log(1e-4), log(0.5), length.out = 60))),
                            n_inits = 5,
                            n_grid_perturb = 3,
                            nstart = 10,
                            alpha_upper = 8.5,
                            selector = "H_elbow",
                            aic_penalty_mult = 1.0,
                            bic_penalty_mult = 1.0,
                            h_elbow_jump_frac = 0.03,
                            h_elbow_jump_min = 0.0,
                            verbose = TRUE) {

  ## Pre-define the grid-perturbation shift factors (safe for any n_grid_perturb)
  shift_factors <- c(0.7, 1.0, 1.4)
  n_grid_perturb <- min(n_grid_perturb, length(shift_factors))   # guard against p > 3

  total_trials <- n_inits + n_grid_perturb
  results <- list(); trial <- 0L

  ## --- Phase 1: seed perturbation trials --------------------------
  for (s in seq_len(n_inits)) {
    trial <- trial + 1L
    if (verbose) cat(sprintf("\nStability trial %d/%d: seed perturbation\n", trial, total_trials))
    res <- tryCatch(
      fit_pmhd_aepd_galaxies(
        x, K_max = K_max, lambda_grid = base_lambda_grid,
        nstart = nstart, alpha_upper = alpha_upper,
        selector = selector,
        aic_penalty_mult = aic_penalty_mult,
        bic_penalty_mult = bic_penalty_mult,
        h_elbow_jump_frac = h_elbow_jump_frac,
        h_elbow_jump_min = h_elbow_jump_min,
        seed = 20250420 + 1000L * s, verbose = FALSE, draw_paths = FALSE),
      error = function(e) {
        if (verbose) cat(sprintf("  [WARNING] seed trial %d failed: %s\n", s, conditionMessage(e)))
        NULL
      })
    if (!is.null(res)) {
      results[[trial]] <- list(
        type = sprintf("seed_%d", s),
        Khat = res$best_K,
        mu = sort(res$best_fit$mu),
        H_path = sqrt(res$best_fit$H2),
        H_final = sqrt(res$best_fit$H2_final %||% res$best_fit$H2),
        lambda = res$best_lambda)
    }
  }

  ## --- Phase 2: lambda-grid perturbation trials -------------------
  for (p in seq_len(n_grid_perturb)) {
    trial <- trial + 1L
    new_max <- 0.5 * shift_factors[p]
    grid_p <- c(0, exp(seq(log(1e-4), log(new_max), length.out = length(base_lambda_grid) - 1L)))
    if (verbose) cat(sprintf("\nStability trial %d/%d: lambda grid up to %.3f\n", trial, total_trials, new_max))
    res <- tryCatch(
      fit_pmhd_aepd_galaxies(
        x, K_max = K_max, lambda_grid = grid_p,
        nstart = nstart, alpha_upper = alpha_upper,
        selector = selector,
        aic_penalty_mult = aic_penalty_mult,
        bic_penalty_mult = bic_penalty_mult,
        h_elbow_jump_frac = h_elbow_jump_frac,
        h_elbow_jump_min = h_elbow_jump_min,
        seed = 20250420 + 10000L * p, verbose = FALSE, draw_paths = FALSE),
      error = function(e) {
        if (verbose) cat(sprintf("  [WARNING] grid-perturbation trial %d failed: %s\n", p, conditionMessage(e)))
        NULL
      })
    if (!is.null(res)) {
      results[[trial]] <- list(
        type = sprintf("grid_pert_%d", p),
        Khat = res$best_K,
        mu = sort(res$best_fit$mu),
        H_path = sqrt(res$best_fit$H2),
        H_final = sqrt(res$best_fit$H2_final %||% res$best_fit$H2),
        lambda = res$best_lambda)
    }
  }

  ## Remove any NULL slots left by failed trials
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) stop("All stability trials failed.")
  results
}

summarize_stability <- function(stab) {
  df <- data.frame(
    trial = sapply(stab, `[[`, "type"),
    Khat = sapply(stab, `[[`, "Khat"),
    H_path = sapply(stab, `[[`, "H_path"),
    H_final = sapply(stab, `[[`, "H_final"),
    lambda = sapply(stab, `[[`, "lambda"),
    stringsAsFactors = FALSE)
  mus_list <- lapply(stab, `[[`, "mu")
  all_equal_K <- length(unique(sapply(mus_list, length))) == 1L
  mu_range <- NA_real_
  if (all_equal_K) {
    mus_mat <- do.call(rbind, mus_list)
    mu_range <- max(apply(mus_mat, 2, function(u) max(u) - min(u)))
  }
  list(table = df,
       all_same_Khat = length(unique(df$Khat)) == 1L,
       mu_range_max = mu_range)
}

# ============================================================

# ---- PMHD-MCP-R refit functions ------------------------------
# ============================================================
# Smoothed-model MHD refit on the asinh scale (Basu & Lindsay, 1994)
# ============================================================
#
# Purpose
# -------
# Post-selection bias-corrected refit for PMHD-MCP. After the
# regularization path has selected (K*, lambda*), refit the K*-component
# AEPD mixture by minimizing the Hellinger distance between
#
#     ghat   (the asinh-scale KDE used throughout the pipeline)
# and
#     A f_Theta = sum_k pi_k * (A g_{theta_k}),
#
# where A is EXACTLY the same linear smoothing operator that produced
# ghat from the data:
#
#     (A f)(x) = [ K_bw  *  f_Y ]( asinh(x/s) ) * dy/dx,
#     f_Y(y)   = f( s*sinh(y) ) * s*cosh(y),     dy/dx = 1/sqrt(x^2+s^2),
#
# optionally followed by the same tail floor (1-eps)*(.) + eps*ref used
# inside estimate_empirical_density_grid(). Because convolution is
# linear and sum_k pi_k = 1, the mixture structure is preserved
# component-wise, so the same MM decomposition as pmhd_mcp_fit_one()
# applies verbatim with lambda = 0 (no MCP penalty, K fixed).
#
# At the model, E[ghat] ~ A g_0, so the criterion is Fisher consistent
# for any fixed bandwidth: the O(bw^2) smoothing bias that inflates the
# fitted scale parameters under the plug-in criterion cancels.
#
# Dependencies (defined in simulation_stability_elbow_parallel_main.R):
#   daepd, trapz_weights, normalize_grid_density, H2_grid,
#   estimate_empirical_density_grid
# Optional: nloptr (BOBYQA); falls back to stats::optim(L-BFGS-B).
#
# Exported functions:
#   make_asinh_smoother(x, kde, pad_bw = 8)
#   smooth_aepd_component(sm, mu, sigma, alpha, tau)
#   smoothed_mixture_grid(sm, pi, mu, sigma, alpha, tau)
#   mhd_smoothed_refit(x, init, kde = NULL, sm = NULL, ...)
# ============================================================

# ------------------------------------------------------------
# Build the smoothing operator A matching the asinh-scale KDE
# ------------------------------------------------------------
# x   : the data vector used to build the KDE (needed only to replicate
#       the t4 tail-floor reference density exactly; pass the same x).
# kde : the object returned by
#       estimate_empirical_density_grid(x, method = "asinh", ...).
# pad_bw : padding of the y-grid, in units of the bandwidth. The model
#       is evaluated analytically on the padded grid, so mass flowing
#       across the data-driven window edges is handled correctly.
make_kde_smoother <- function(x, kde, pad_bw = 8, tail_eps = NULL,
                              ref_df = 4) {
  std <- !is.null(kde$method) && kde$method == "standard"
  s  <- if (std) NA_real_ else kde$scale
  bw <- kde$bw
  grid <- kde$grid
  G  <- length(grid)

  # Reconstruct the uniform working grid used by density():
  #   asinh method: y = asinh(x/s);  standard method: y = x.
  y  <- if (std) grid else asinh(grid / s)
  dy <- diff(y)
  if (max(dy) - min(dy) > 1e-8 * mean(dy))
    stop("make_kde_smoother: working grid is not uniform; ",
         "kde does not look like the output of ",
         "estimate_empirical_density_grid().")
  dy <- mean(dy)

  # Padded y-grid (model evaluated analytically there).
  npad  <- as.integer(ceiling(pad_bw * bw / dy))
  M     <- npad                       # kernel half-width in grid steps
  y_pad <- seq(y[1L] - npad * dy, by = dy, length.out = G + 2L * npad)
  Gp    <- length(y_pad)
  if (std) {
    x_pad <- y_pad
    dxdy_pad <- rep(1, Gp)
  } else {
    x_pad <- s * sinh(y_pad)
    dxdy_pad <- s * cosh(y_pad)       # f_Y(y) = f_X(x(y)) * dx/dy
  }

  # Discrete Gaussian kernel, normalized so that sum(kk)*dy = 1
  # (matches the Gaussian kernel of density() up to discretization,
  #  with the truncation at +-pad_bw*bw corrected by renormalization).
  kk <- dnorm(seq.int(-M, M) * dy, sd = bw)
  kk <- kk / (sum(kk) * dy)

  # Banded smoothing matrix S (G x Gp):
  #   (S %*% f_pad)[i] = sum_m kk[(i + npad - m) + M + 1] * f_pad[m] * dy
  # i.e. the convolution evaluated at the i-th ORIGINAL grid point.
  ii <- matrix(seq_len(G),  nrow = G, ncol = Gp)
  mm <- matrix(seq_len(Gp), nrow = G, ncol = Gp, byrow = TRUE)
  dd <- ii + npad - mm
  S  <- matrix(0, G, Gp)
  ok <- abs(dd) <= M
  S[ok] <- kk[dd[ok] + M + 1L] * dy
  rm(ii, mm, dd, ok)

  jac_x <- if (std) rep(1, G) else 1 / sqrt(grid^2 + s^2)  # dy/dx

  # Replicate the tail floor of estimate_empirical_density_grid():
  #   ghat <- (1-eps)*ghat + eps*ref,  ref = t_{ref_df}((x-loc)/sc)/sc.
  # Since sum_k pi_k = 1, applying the same affine map per component
  # keeps the mixture structure intact.
  if (is.null(tail_eps)) tail_eps <- 1e-6   # pipeline default
  eps <- 0; ref_g <- NULL
  if (is.finite(tail_eps) && tail_eps > 0) {
    loc <- median(x)
    sc  <- mad(x, constant = 1.4826)
    if (!is.finite(sc) || sc <= 0) sc <- sd(x)
    if (!is.finite(sc) || sc <= 0) sc <- s
    ref_g <- dt((grid - loc) / sc, df = ref_df) / sc
    eps   <- min(max(tail_eps, 0), 0.05)
  }

  list(grid = grid, wq = kde$wq, ghat = kde$ghat,
       s = s, bw = bw, dy = dy, npad = npad, M = M,
       x_pad = x_pad, dxdy_pad = dxdy_pad,
       S = S, jac_x = jac_x, eps = eps, ref_g = ref_g)
}

# ------------------------------------------------------------
# A g_theta : smoothed AEPD component on the original x-grid
# ------------------------------------------------------------
smooth_aepd_component <- function(sm, mu, sigma, alpha, tau) {
  fy <- daepd(sm$x_pad, mu, sigma, alpha, tau) * sm$dxdy_pad
  fs <- drop(sm$S %*% fy) * sm$jac_x
  if (sm$eps > 0) fs <- (1 - sm$eps) * fs + sm$eps * sm$ref_g
  pmax(fs, 1e-300)
}

smoothed_mixture_grid <- function(sm, pi, mu, sigma, alpha, tau) {
  f <- rep(0, length(sm$grid))
  for (k in seq_along(pi))
    f <- f + pi[k] * smooth_aepd_component(sm, mu[k], sigma[k],
                                           alpha[k], tau[k])
  pmax(f, 1e-300)
}

# ------------------------------------------------------------
# Box-constrained local optimizer: nloptr::BOBYQA if available,
# otherwise stats::optim(L-BFGS-B). Same reparametrization as
# pmhd_mcp_fit_one: (mu, log sigma, log alpha, logit tau).
# ------------------------------------------------------------
.refit_optimize <- function(par0, fn, lb, ub, maxeval = 300) {
  if (requireNamespace("nloptr", quietly = TRUE)) {
    res <- tryCatch(
      nloptr::nloptr(par0, fn, lb = lb, ub = ub,
                     opts = list(algorithm = "NLOPT_LN_BOBYQA",
                                 xtol_rel = 1e-6, maxeval = maxeval)),
      error = function(e) NULL)
    if (!is.null(res))
      return(list(par = res$solution, value = res$objective))
  }
  par0 <- pmin(pmax(par0, lb + 1e-8), ub - 1e-8)
  res <- tryCatch(
    stats::optim(par0, fn, method = "L-BFGS-B",
                 lower = lb, upper = ub,
                 control = list(maxit = maxeval)),
    error = function(e) NULL)
  if (!is.null(res)) return(list(par = res$par, value = res$value))
  list(par = par0, value = fn(par0))
}

# ------------------------------------------------------------
# Smoothed-model MHD refit at fixed K (lambda = 0), ANCHORED
# ------------------------------------------------------------
# Criterion:  J(Theta) = H2(ghat, A f_Theta) + gamma * H2(ghat, f_Theta)
#
# The smoothed term removes the plug-in smoothing bias (Fisher
# consistent at fixed bandwidth); the small raw-criterion anchor
# (gamma_anchor, default 0.10) restores the repulsion of degenerate
# spike components: under pure model smoothing, a collapsed component
# (sigma -> 0) convolves to a kernel-width bump and is barely
# penalized, so with noisy ghat / poor warm starts the criterion is
# nearly flat for sigma below the local effective bandwidth. The
# anchor re-creates the cliff (raw affinity of a spike -> 0) at the
# cost of retaining only a fraction gamma/(1+gamma) of the original
# bias -- i.e. ~90% of the bias correction is kept at gamma = 0.10.
#
# Additional safeguards:
#   * per-component lower bound sigma_k >= shrink_floor * sigma_init_k
#     (refit may shrink scales -- that is the point -- but in these
#     designs the true/raw-init scale ratio is ~0.6-0.75, so allowing
#     up to 75% shrink is generous while blocking spike collapse);
#   * acceptance guard: the refit is returned only if it improves the
#     anchored objective over the warm start, otherwise the warm start
#     is returned unchanged (accepted = FALSE).
mhd_smoothed_refit <- function(x, init, kde = NULL, sm = NULL,
                               bw_adjust = 1.0, grid_n = 512,
                               density_method = "asinh",
                               density_tail_eps = 1e-6,
                               tol = 1e-6, max_iter = 100,
                               maxeval = 300,
                               alpha_lower = 1.0, alpha_upper = 3.0,
                               tau_lower = 0.01, tau_upper = 0.99,
                               gamma_anchor = 0.10,
                               shrink_floor = 0.25,
                               raw_trust = 0,
                               accept_guard = TRUE,
                               smooth_model = TRUE) {
  if (is.null(kde))
    kde <- estimate_empirical_density_grid(
      x, grid_n = grid_n, bw_adjust = bw_adjust,
      method = density_method, tail_eps = density_tail_eps)
  if (is.null(sm) && smooth_model)
    sm <- make_kde_smoother(x, kde, tail_eps = density_tail_eps)

  grid <- kde$grid; wq <- kde$wq; ghat <- kde$ghat
  if (!smooth_model) gamma_anchor <- 0
  gam <- max(gamma_anchor, 0)

  comp_s <- if (smooth_model) {
    function(mu, sg, al, tt) smooth_aepd_component(sm, mu, sg, al, tt)
  } else {
    function(mu, sg, al, tt) daepd(grid, mu, sg, al, tt)
  }
  comp_r <- function(mu, sg, al, tt) daepd(grid, mu, sg, al, tt)

  K <- length(init$pi)
  pi_k <- pmax(init$pi, 1e-8); pi_k <- pi_k / sum(pi_k)
  mu_k <- init$mu; sigma_k <- init$sigma
  alpha_k <- pmin(pmax(init$alpha, alpha_lower), alpha_upper)
  tau_k   <- pmin(pmax(init$tau, tau_lower), tau_upper)

  sigma_glob_lo <- 1e-3 * sd(x); sigma_up <- 3 * diff(range(x))
  # Median-based absolute floor: shrink_floor * sigma_init is useless if
  # the init itself carries a near-degenerate (tiny-pi, tiny-sigma)
  # component that survived path trimming -- the pi ~ B^2 update can
  # then AMPLIFY the spike. Floor every scale at 10% of the median init
  # scale, and repair clearly degenerate init components by resetting
  # them to the median scale so the MM responsibilities can reallocate.
  sig_med <- stats::median(pmax(init$sigma, 0), na.rm = TRUE)
  if (!is.finite(sig_med) || sig_med <= 0) sig_med <- sd(x) / 4
  degk <- which(sigma_k < 0.05 * sig_med)
  if (length(degk)) sigma_k[degk] <- sig_med
  sigma_lo_k <- pmax(sigma_glob_lo,
                     shrink_floor * pmax(init$sigma, 0),
                     0.10 * sig_med)
  sigma_k <- pmin(pmax(sigma_k, sigma_lo_k * 1.001), sigma_up / 1.001)

  Fs <- matrix(0, length(grid), K)   # smoothed components
  Fr <- matrix(0, length(grid), K)   # raw components
  set_comp <- function(k) {
    Fs[, k] <<- comp_s(mu_k[k], sigma_k[k], alpha_k[k], tau_k[k])
    Fr[, k] <<- comp_r(mu_k[k], sigma_k[k], alpha_k[k], tau_k[k])
  }
  for (k in seq_len(K)) set_comp(k)

  # anchored affinity (larger is better); J = (1+gam) - A_total
  A_total <- function() {
    fs <- pmax(drop(Fs %*% pi_k), 1e-300)
    fr <- pmax(drop(Fr %*% pi_k), 1e-300)
    sum(sqrt(ghat * fs) * wq) + gam * sum(sqrt(ghat * fr) * wq)
  }
  A_init <- A_total(); A_old <- A_init

  # raw-criterion trust region: a legitimate bias correction increases
  # the RAW Hellinger distance to ghat only slightly (the raw optimum is
  # the warm start); a spiked/misplaced component increases it a lot.
  raw_H2_now <- function()
    H2_grid(ghat, normalize_grid_density(
      pmax(drop(Fr %*% pi_k), 1e-300), wq), wq)
  rawH2_init <- raw_H2_now()

  pi0_s <- pi_k; mu0_s <- mu_k; sg0_s <- sigma_k
  al0_s <- alpha_k; ta0_s <- tau_k     # warm-start snapshot

  iter_used <- 0L
  for (iter in seq_len(max_iter)) {
    iter_used <- iter
    fs_mix <- pmax(drop(Fs %*% pi_k), 1e-300)
    fr_mix <- pmax(drop(Fr %*% pi_k), 1e-300)
    Rs <- sweep(sweep(Fs, 2, pi_k, "*"), 1, fs_mix, "/")
    Rr <- sweep(sweep(Fr, 2, pi_k, "*"), 1, fr_mix, "/")

    for (k in seq_len(K)) {
      gs_k <- ghat * Rs[, k]
      gr_k <- ghat * Rr[, k]
      obj_theta <- function(par) {
        mu <- par[1]; sg <- exp(par[2]); al <- exp(par[3]); tt <- plogis(par[4])
        fs <- comp_s(mu, sg, al, tt)
        a  <- sum(sqrt(gs_k * fs) * wq)
        if (gam > 0) {
          fr <- comp_r(mu, sg, al, tt)
          a <- a + gam * sum(sqrt(gr_k * fr) * wq)
        }
        -a
      }
      lb <- c(-Inf, log(sigma_lo_k[k]), log(alpha_lower), qlogis(tau_lower))
      ub <- c( Inf, log(sigma_up),      log(alpha_upper), qlogis(tau_upper))
      res <- .refit_optimize(
        c(mu_k[k], log(sigma_k[k]), log(alpha_k[k]), qlogis(tau_k[k])),
        obj_theta, lb, ub, maxeval = maxeval)
      mu_k[k]    <- res$par[1]
      sigma_k[k] <- exp(res$par[2])
      alpha_k[k] <- exp(res$par[3])
      tau_k[k]   <- plogis(res$par[4])
      set_comp(k)
    }

    # lambda = 0: pi_k proportional to B_k^2 (anchored B)
    B_k <- vapply(seq_len(K), function(k)
      2 * (sum(sqrt(pmax(ghat * Rs[, k] * Fs[, k], 0)) * wq) +
           gam * sum(sqrt(pmax(ghat * Rr[, k] * Fr[, k], 0)) * wq)),
      numeric(1))
    pi_k <- B_k^2
    if (sum(pi_k) <= 0) pi_k <- rep(1 / K, K) else pi_k <- pi_k / sum(pi_k)

    A_new <- A_total()
    if (abs(A_new - A_old) < tol) break
    A_old <- A_new
  }

  A_final <- A_total()
  rawH2_final <- raw_H2_now()
  accepted <- TRUE
  raw_tol <- max(raw_trust, 1e-10)
  trust_ok <- is.finite(rawH2_final) &&
    (rawH2_final - rawH2_init) <= raw_tol
  if (isTRUE(accept_guard) &&
      !(is.finite(A_final) && A_final >= A_init && trust_ok)) {
    pi_k <- pi0_s; mu_k <- mu0_s; sigma_k <- sg0_s
    alpha_k <- al0_s; tau_k <- ta0_s
    rawH2_final <- rawH2_init
    accepted <- FALSE
  }

  fsm <- if (smooth_model)
    normalize_grid_density(
      smoothed_mixture_grid(sm, pi_k, mu_k, sigma_k, alpha_k, tau_k), wq)
  else {
    fm <- rep(0, length(grid))
    for (k in seq_len(K)) fm <- fm + pi_k[k] * Fr[, k]
    normalize_grid_density(pmax(fm, 1e-300), wq)
  }

  ord <- order(mu_k)
  list(pi = pi_k[ord], mu = mu_k[ord], sigma = sigma_k[ord],
       alpha = alpha_k[ord], tau = tau_k[ord],
       Khat = K, H2_smoothed = H2_grid(ghat, fsm, wq),
       iter = iter_used, accepted = accepted,
       raw_H2_increase = rawH2_final - rawH2_init,
       gamma_anchor = gam, smoothed = smooth_model)
}

# ============================================================
# Integration into scenario1.R (replace the post-cache selection)
# ============================================================
# select_pmhd_from_cache <- function(pm, x = NULL, kde = NULL, sm = NULL) {
#   ...                                    # unchanged up to `fit`
#   if (is.null(fit)) return(op)
#   if (!is.null(x) && is.finite(fit$Khat) && fit$Khat >= 1L) {
#     refit <- tryCatch(mhd_smoothed_refit(
#       x, init = fit[c("pi","mu","sigma","alpha","tau")],
#       kde = kde, sm = sm,
#       alpha_lower = alpha_lower, alpha_upper = alpha_upper),
#       error = function(e) NULL)
#     if (!is.null(refit)) fit <- refit
#   }
#   fh <- normalize_grid_density(aepd_density_grid(grid, fit), wq)
#   ...                                    # metrics unchanged: fh is the
# }                                        # RAW AEPD mixture at refit params
#
# br <- lapply(seq_along(br), function(b) {
#   r <- br[[b]]
#   set.seed(seed + b + 1000L * ni)        # same stream as cache generation
#   x_b <- rmix_aepd(n, pi0, mu0, sigma0, alpha0, tau0)
#   kde_b <- estimate_empirical_density_grid(
#     x_b, grid_n = 512, bw_adjust = bw_adjust,
#     method = density_method, tail_eps = density_tail_eps)
#   sm_b <- make_asinh_smoother(x_b, kde_b, tail_eps = density_tail_eps)
#   r$pmhd <- select_pmhd_from_cache(r$pmhd_cache, x = x_b,
#                                    kde = kde_b, sm = sm_b)
#   r
# })
#
# Notes
# -----
# 1. The metrics (H, ISE, MSE) must be computed from the RAW AEPD
#    mixture at the refit parameters (aepd_density_grid), exactly as
#    before. The smoothing lives only inside the refit criterion.
# 2. The KDE settings passed to estimate_empirical_density_grid must be
#    identical to those used inside select_model_stability (bw_adjust,
#    grid_n, method, tail_eps), so that ghat and the smoother are the
#    same operator bit-for-bit.
# 3. Selection (P(Khat = K0), Kbar) is untouched: Khat comes from the
#    path; the refit holds K fixed.
# 4. For method = "standard" the same machinery applies with the
#    identity transform: y = x, dxdy = 1, jac_x = 1, and the uniform
#    x-grid of density(); only make_asinh_smoother's transform lines
#    change.

# Back-compatibility wrapper (asinh-only name used in earlier scripts)
make_asinh_smoother <- function(x, kde, ...) make_kde_smoother(x, kde, ...)

# ---- Final galaxies workflow ---------------------------------
write_pmhd_path_outputs <- function(pmhd, output_dir, x) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  write.csv(pmhd$path_df,
            file.path(output_dir, "galaxies_selector_path.csv"),
            row.names = FALSE)
  write.csv(pmhd$compressed_df,
            file.path(output_dir, "galaxies_compressed_path.csv"),
            row.names = FALSE)
  selector_summary <- data.frame(
    method = c("H_elbow", "H-AIC", "H-BIC"),
    Khat = c(pmhd$K_H_elbow, pmhd$K_aic, pmhd$K_bic),
    lambda = c(pmhd$lambda_H_elbow, pmhd$lambda_aic, pmhd$lambda_bic),
    stringsAsFactors = FALSE)
  write.csv(selector_summary,
            file.path(output_dir, "galaxies_selector_summary.csv"),
            row.names = FALSE)
  plot_selector_paths(pmhd,
                      file = file.path(output_dir, "galaxies_selector_paths.pdf"))
  saveRDS(list(pmhd = pmhd,
               data_summary = list(n = length(x), range = range(x)),
               selector_summary = selector_summary),
          file.path(output_dir, "galaxies_results.rds"))
  invisible(selector_summary)
}

has_galaxies_final_outputs <- function(dir) {
  required <- c("galaxies_pmhd_refit_summary.csv",
                "galaxies_pmhd_refit_table.csv",
                "galaxies_pmhd_refit_fit.pdf",
                "galaxies_pmhd_refit_results.rds",
                "galaxies_selector_paths.pdf")
  dir.exists(dir) && all(file.exists(file.path(dir, required)))
}

sync_galaxies_final_outputs <- function(
    output_dir,
    source_dirs = c(file.path(galaxies_dir, "K6final"),
                    file.path(galaxies_dir, "K6"))) {
  source_dirs <- normalizePath(source_dirs, winslash = "/", mustWork = FALSE)
  source_dir <- source_dirs[vapply(source_dirs, has_galaxies_final_outputs,
                                   logical(1L))][1L]
  if (is.na(source_dir)) {
    stop("No complete source results found. Expected one of: ",
         paste(source_dirs, collapse = ", "))
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  files <- list.files(source_dir, full.names = FALSE)
  for (f in files) {
    file.copy(file.path(source_dir, f), file.path(output_dir, f),
              overwrite = TRUE)
  }
  cat("Copied final galaxies outputs from:\n  ", source_dir,
      "\nto:\n  ", normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      "\n", sep = "")
  summary_file <- file.path(output_dir, "galaxies_pmhd_refit_summary.csv")
  if (file.exists(summary_file)) print(read.csv(summary_file))
  invisible(readRDS(file.path(output_dir, "galaxies_pmhd_refit_results.rds")))
}

galaxies_paper_lambda <- 0.0156424307087534

galaxies_paper_pmhd_init <- function() {
  list(
    pi = c(0.088655, 0.295435, 0.188174,
           0.237677, 0.152450, 0.037610),
    mu = c(9667.573161, 19442.665510, 19637.943939,
           20556.806756, 23220.191817, 32711.028313),
    sigma = c(649.475270, 919.743497, 769.579206,
              1122.924641, 714.683341, 939.821831),
    alpha = c(2.063700, 1.993899, 1.260172,
              2.430273, 2.418036, 2.790063),
    tau = c(0.485607, 0.273144, 0.431243,
            0.258206, 0.385829, 0.412058)
  )
}

fit_pmhd_fixed_paper <- function(x,
                                 lambda = galaxies_paper_lambda,
                                 gamma = 3,
                                 grid_n = 1024,
                                 bw_adjust = 0.80,
                                 tol = 1e-6,
                                 max_iter = 80,
                                 maxeval = 200,
                                 alpha_lower = 1.0,
                                 alpha_upper = 8.5,
                                 seed = 20250420,
                                 delta_final = NULL,
                                 bic_params_per_K = 4L) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(delta_final)) delta_final <- max(0.001, 2 / length(x))

  bw0 <- bw.nrd0(x) * bw_adjust
  dens <- density(x, n = grid_n, bw = bw0,
                  from = min(x) - 5 * bw0, to = max(x) + 5 * bw0)
  grid <- dens$x
  wq <- trapz_weights(grid)
  ghat <- normalize_grid_density(dens$y, wq)
  init <- galaxies_paper_pmhd_init()

  cat(sprintf(
    "Fast mode: fixed K=6, lambda=%.8f, paper PMHD-MCP warm start\n",
    lambda))
  fit <- pmhd_mcp_fit_one(
    x, K_init = length(init$pi), lambda = lambda, gamma = gamma,
    grid = grid, wq = wq, ghat = ghat,
    tol = tol, max_iter = max_iter, delta_final = delta_final,
    init = init, maxeval = maxeval,
    alpha_lower = alpha_lower, alpha_upper = alpha_upper
  )
  fit <- sort_fit_by_mu(fit)

  Khat <- fit$Khat
  H2 <- fit$H2
  H <- sqrt(H2)
  df <- Khat * bic_params_per_K + (Khat - 1)
  path_df <- data.frame(
    idx = 1L,
    lambda = lambda,
    log10_lambda = log10(lambda),
    K_raw = Khat,
    K_mono = Khat,
    H2 = H2,
    H = H,
    nH2 = length(x) * H2,
    df = df,
    score_aic = length(x) * H2 + 2 * df,
    score_bic = length(x) * H2 + log(length(x)) * df
  )
  compressed_df <- compress_lambda_path(lambda, Khat, H2, H)
  h_elbow <- choose_H_elbow(compressed_df, n_obs = length(x))

  list(
    best_fit = fit,
    best_lambda = lambda,
    best_K = Khat,
    best_idx = 1L,
    selector = "fixed_paper_lambda",
    path_df = path_df,
    compressed_df = compressed_df,
    fits_by_lambda = list(fit),
    fit_H_elbow = fit,
    K_H_elbow = Khat,
    lambda_H_elbow = lambda,
    idx_H_elbow = 1L,
    H_elbow_details = h_elbow,
    fit_aic = fit,
    K_aic = Khat,
    lambda_aic = lambda,
    idx_aic = 1L,
    fit_bic = fit,
    K_bic = Khat,
    lambda_bic = lambda,
    idx_bic = 1L,
    kde = list(method = "standard", grid = grid, wq = wq,
               ghat = ghat, bw = bw0),
    config = list(mode = "fixed_paper_warm_start",
                  lambda = lambda, K_fixed = length(init$pi),
                  bw_adjust = bw_adjust, grid_n = grid_n,
                  max_iter = max_iter, maxeval = maxeval,
                  delta_final = delta_final)
  )
}

## ---- KDE helper expected by smoothed_refit.R ------------------
estimate_empirical_density_grid <- function(x, grid_n = 1024,
                                            bw_adjust = 0.80,
                                            method = "standard",
                                            tail_eps = 0,
                                            from = NULL, to = NULL) {
  method <- match.arg(method, c("standard"))
  bw0 <- bw.nrd0(x) * bw_adjust
  if (is.null(from)) from <- min(x) - 5 * bw0
  if (is.null(to))   to   <- max(x) + 5 * bw0
  dens <- density(x, n = grid_n, bw = bw0, from = from, to = to)
  grid <- dens$x
  wq <- trapz_weights(grid)
  ghat <- normalize_grid_density(dens$y, wq)
  if (is.finite(tail_eps) && tail_eps > 0) {
    loc <- median(x)
    sc <- mad(x, constant = 1.4826)
    if (!is.finite(sc) || sc <= 0) sc <- sd(x)
    ref <- dt((grid - loc) / sc, df = 4) / sc
    eps <- min(max(tail_eps, 0), 0.05)
    ghat <- normalize_grid_density((1 - eps) * ghat + eps * ref, wq)
  }
  list(method = method, grid = grid, wq = wq, ghat = ghat, bw = bw0)
}

sort_fit_by_mu <- function(fit) {
  if (is.null(fit$Khat)) fit$Khat <- length(fit$pi)
  ord <- order(fit$mu)
  fit$pi <- fit$pi[ord]
  fit$mu <- fit$mu[ord]
  fit$sigma <- fit$sigma[ord]
  fit$alpha <- fit$alpha[ord]
  fit$tau <- fit$tau[ord]
  fit
}

sort_gm_by_mu <- function(fit) {
  ord <- order(fit$mu)
  fit$pi <- fit$pi[ord]
  fit$mu <- fit$mu[ord]
  fit$sigma <- fit$sigma[ord]
  fit
}

fit_to_table <- function(fit, label) {
  fit <- sort_fit_by_mu(fit)
  data.frame(
    method = label,
    k = seq_len(fit$Khat),
    pi_hat = round(fit$pi, 6),
    mu_hat = round(fit$mu, 6),
    sigma_hat = round(fit$sigma, 6),
    alpha_hat = round(fit$alpha, 6),
    tau_hat = round(fit$tau, 6),
    shape = sapply(seq_len(fit$Khat), function(k)
      describe_shape(fit$alpha[k], fit$tau[k])),
    stringsAsFactors = FALSE
  )
}

raw_mixture_on_grid <- function(grid, fit) {
  eval_aepd_mixture(grid, fit$pi, fit$mu, fit$sigma,
                    fit$alpha, fit$tau)
}

find_default_results_file <- function(output_dir) {
  file.path(output_dir, "galaxies_results.rds")
}

wait_for_file <- function(file, wait_seconds = 0, poll_seconds = 10) {
  if (file.exists(file)) return(TRUE)
  wait_seconds <- max(0, wait_seconds)
  if (wait_seconds <= 0) return(FALSE)

  deadline <- Sys.time() + wait_seconds
  cat("Waiting for new results file: ", file, "\n", sep = "")
  while (Sys.time() < deadline) {
    Sys.sleep(min(poll_seconds, as.numeric(deadline - Sys.time(),
                                           units = "secs")))
    if (file.exists(file)) return(TRUE)
  }
  FALSE
}

can_reselect_pmhd_path <- function(pmhd) {
  !is.null(pmhd) &&
    !is.null(pmhd$path_df) &&
    !is.null(pmhd$fits_by_lambda) &&
    !is.null(pmhd$kde) &&
    !is.null(pmhd$kde$grid) &&
    !is.null(pmhd$best_fit) &&
    !is.null(pmhd$best_fit$H2) &&
    all(c("lambda", "K_mono", "H2", "H") %in% names(pmhd$path_df))
}

fit_khat <- function(fit) {
  if (is.null(fit) || is.null(fit$Khat)) NA_integer_ else fit$Khat
}

reselect_pmhd_from_path <- function(pmhd, n_obs,
                                    h_elbow_jump_frac = 0.03,
                                    h_elbow_jump_min = 0.0) {
  if (!can_reselect_pmhd_path(pmhd))
    stop("PMHD result does not contain path_df/fits_by_lambda for reselection.")

  path <- pmhd$path_df
  lambda_grid <- path$lambda
  compressed_df <- compress_lambda_path(lambda_grid, path$K_mono,
                                        path$H2, path$H)
  he <- choose_H_elbow(compressed_df, n_obs = n_obs,
                       jump_frac = h_elbow_jump_frac,
                       jump_min = h_elbow_jump_min)
  idx_H_elbow <- he$best_idx
  idx_aic <- if ("score_aic" %in% names(path)) which.min(path$score_aic) else idx_H_elbow
  idx_bic <- if ("score_bic" %in% names(path)) which.min(path$score_bic) else idx_H_elbow
  best_idx <- idx_H_elbow

  pmhd$best_idx <- best_idx
  pmhd$best_fit <- pmhd$fits_by_lambda[[best_idx]]
  pmhd$best_lambda <- lambda_grid[best_idx]
  pmhd$best_K <- fit_khat(pmhd$best_fit)
  pmhd$selector <- "H_elbow"

  pmhd$fit_H_elbow <- pmhd$fits_by_lambda[[idx_H_elbow]]
  pmhd$K_H_elbow <- fit_khat(pmhd$fit_H_elbow)
  pmhd$lambda_H_elbow <- lambda_grid[idx_H_elbow]
  pmhd$idx_H_elbow <- idx_H_elbow
  pmhd$H_elbow_details <- he

  pmhd$fit_aic <- pmhd$fits_by_lambda[[idx_aic]]
  pmhd$K_aic <- fit_khat(pmhd$fit_aic)
  pmhd$lambda_aic <- lambda_grid[idx_aic]
  pmhd$idx_aic <- idx_aic

  pmhd$fit_bic <- pmhd$fits_by_lambda[[idx_bic]]
  pmhd$K_bic <- fit_khat(pmhd$fit_bic)
  pmhd$lambda_bic <- lambda_grid[idx_bic]
  pmhd$idx_bic <- idx_bic

  pmhd$compressed_df <- compressed_df
  pmhd$config$h_elbow_jump_frac <- h_elbow_jump_frac
  pmhd$config$h_elbow_jump_min <- h_elbow_jump_min
  pmhd
}

get_or_fit_pmhd <- function(x, results_file,
                            K_max = 8,
                            lambda_grid = c(0, exp(seq(log(1e-4), log(0.5),
                                                     length.out = 60))),
                            gamma = 3,
                            nstart = 20,
                            grid_n = 1024,
                            bw_adjust = 0.80,
                            tol = 1e-6,
                            max_iter = 200,
                            maxeval = 400,
                            alpha_lower = 1.0,
                            alpha_upper = 8.5,
                            h_elbow_jump_frac = 0.03,
                            h_elbow_jump_min = 0.0,
                            seed = 20250420,
                            force_path_refit = TRUE) {
  if (isTRUE(force_path_refit)) {
    return(fit_pmhd_aepd_galaxies(
      x,
      K_max = K_max,
      lambda_grid = lambda_grid,
      gamma = gamma,
      nstart = nstart,
      grid_n = grid_n,
      bw_adjust = bw_adjust,
      tol = tol,
      max_iter = max_iter,
      maxeval = maxeval,
      alpha_lower = alpha_lower,
      alpha_upper = alpha_upper,
      selector = "H_elbow",
      h_elbow_jump_frac = h_elbow_jump_frac,
      h_elbow_jump_min = h_elbow_jump_min,
      seed = seed,
      draw_paths = FALSE,
      verbose = TRUE
    ))
  }

  if (!file.exists(results_file)) {
    stop("New PMHD results file not found: ", results_file,
         "\nRun galaxies_analysis_Helbow_framework_manualplot.R first, ",
         "or call run_galaxies_final_analysis(force_path_refit = TRUE).")
  }

  obj <- readRDS(results_file)
  if (!can_reselect_pmhd_path(obj$pmhd)) {
    stop("New PMHD results file does not contain a reusable path: ",
         results_file)
  }
  cat("Using new PMHD path from: ", results_file, "\n", sep = "")
  cat("Reselecting order with good_init jump selector.\n")
  reselect_pmhd_from_path(
    obj$pmhd, n_obs = length(x),
    h_elbow_jump_frac = h_elbow_jump_frac,
    h_elbow_jump_min = h_elbow_jump_min
  )
}

run_galaxies_final_analysis <- function(
    output_dir = Sys.getenv("GALAXIES_FINAL_OUTPUT_DIR", unset = file.path(galaxies_dir, "K6final")),
    results_file = Sys.getenv("GALAXIES_RESULTS_FILE",
                              unset = find_default_results_file(output_dir)),
    wait_for_results = FALSE,
    wait_seconds = 0,
    poll_seconds = 10,
    force_path_refit = TRUE,
    fast_paper_mode = FALSE,
    fixed_lambda = galaxies_paper_lambda,
    use_locked_results = TRUE,
    K_max = 8,
    lambda_n = 61,
    nstart = 20,
    seed = 20250420,
    gamma = 3,
    grid_n = 1024,
    bw_adjust = 0.80,
    alpha_lower = 1.0,
    alpha_upper = 8.5,
    tol = 1e-6,
    max_iter = 200,
    maxeval = 400,
    refit_maxeval = 300,
    refit_max_iter = 100,
    gamma_anchor = 0.10,
    shrink_floor = 0.25,
    raw_trust = 0,
    h_elbow_jump_frac = 0.03,
    h_elbow_jump_min = 0.0) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  data(galaxies, package = "MASS")
  x <- as.numeric(galaxies)
  if (isTRUE(use_locked_results)) {
    return(sync_galaxies_final_outputs(
      output_dir = output_dir,
      source_dirs = c(file.path(galaxies_dir, "K6"),
                      file.path(galaxies_dir, "K6final"))
    ))
  }
  lambda_grid <- c(0, exp(seq(log(1e-4), log(0.5),
                              length.out = lambda_n - 1L)))

  if (!isTRUE(force_path_refit) &&
      isTRUE(wait_for_results) &&
      !wait_for_file(results_file, wait_seconds, poll_seconds)) {
    stop("Timed out waiting for new PMHD results file: ", results_file)
  }

  pmhd <- if (isTRUE(fast_paper_mode)) {
    fit_pmhd_fixed_paper(
      x, lambda = fixed_lambda, gamma = gamma,
      grid_n = grid_n, bw_adjust = bw_adjust,
      tol = tol, max_iter = min(max_iter, 80L),
      maxeval = min(maxeval, 200L),
      alpha_lower = alpha_lower, alpha_upper = alpha_upper,
      seed = seed
    )
  } else {
    get_or_fit_pmhd(
      x, results_file = results_file,
      K_max = K_max, lambda_grid = lambda_grid, gamma = gamma,
      nstart = nstart, grid_n = grid_n, bw_adjust = bw_adjust,
      tol = tol, max_iter = max_iter, maxeval = maxeval,
      alpha_lower = alpha_lower, alpha_upper = alpha_upper,
      h_elbow_jump_frac = h_elbow_jump_frac,
      h_elbow_jump_min = h_elbow_jump_min,
      seed = seed, force_path_refit = force_path_refit
    )
  }

  write_pmhd_path_outputs(pmhd, output_dir, x)

  raw_fit <- sort_fit_by_mu(pmhd$best_fit)

  kde <- estimate_empirical_density_grid(
    x, grid_n = grid_n, bw_adjust = bw_adjust,
    method = "standard", tail_eps = 0,
    from = min(pmhd$kde$grid), to = max(pmhd$kde$grid)
  )
  sm <- make_asinh_smoother(x, kde, tail_eps = 0)

  cat(sprintf("Refitting PMHD at fixed K=%d from lambda=%.6g\n",
              raw_fit$Khat, pmhd$best_lambda))
  t0 <- proc.time()[["elapsed"]]
  refit <- mhd_smoothed_refit(
    x,
    init = raw_fit[c("pi", "mu", "sigma", "alpha", "tau")],
    kde = kde,
    sm = sm,
    tol = tol,
    max_iter = refit_max_iter,
    maxeval = refit_maxeval,
    alpha_lower = alpha_lower,
    alpha_upper = alpha_upper,
    gamma_anchor = gamma_anchor,
    shrink_floor = shrink_floor,
    raw_trust = raw_trust,
    smooth_model = FALSE
  )
  refit_time <- proc.time()[["elapsed"]] - t0
  if (is.null(refit$Khat)) refit$Khat <- length(refit$pi)
  refit <- sort_fit_by_mu(refit)

  raw_on_grid <- normalize_grid_density(
    raw_mixture_on_grid(kde$grid, raw_fit), kde$wq)
  refit_on_grid <- normalize_grid_density(
    raw_mixture_on_grid(kde$grid, refit), kde$wq)
  raw_H_check <- sqrt(H2_grid(kde$ghat, raw_on_grid, kde$wq))
  refit_H <- sqrt(H2_grid(kde$ghat, refit_on_grid, kde$wq))

  gm <- fit_gm_mclust(x, Kmax = K_max)
  gm <- sort_gm_by_mu(gm)
  gm_on_grid <- normalize_grid_density(
    eval_gm(kde$grid, gm$pi, gm$mu, gm$sigma), kde$wq)
  gm_H <- sqrt(H2_grid(kde$ghat, gm_on_grid, kde$wq))

  table_out <- rbind(
    fit_to_table(raw_fit, "PMHD-MCP raw"),
    fit_to_table(refit, "PMHD-MCP-R")
  )
  table_file <- file.path(output_dir, "galaxies_pmhd_refit_table.csv")
  write.csv(table_out, table_file, row.names = FALSE)

  summary_df <- data.frame(
    method = c("PMHD-MCP raw", "PMHD-MCP-R", "Gaussian mixture"),
    Khat = c(raw_fit$Khat, refit$Khat, gm$Khat),
    H = c(raw_H_check, refit_H, gm_H),
    lambda = c(pmhd$best_lambda, NA_real_, NA_real_),
    refit_accepted = c(NA, isTRUE(refit$accepted), NA),
    stringsAsFactors = FALSE
  )
  summary_file <- file.path(output_dir, "galaxies_pmhd_refit_summary.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)

  plot_grid <- seq(min(x) - 2000, max(x) + 2000, length.out = 1000)
  f_raw <- normalize_grid_density(raw_mixture_on_grid(plot_grid, raw_fit),
                                  trapz_weights(plot_grid))
  f_refit <- normalize_grid_density(raw_mixture_on_grid(plot_grid, refit),
                                    trapz_weights(plot_grid))
  f_gm <- normalize_grid_density(eval_gm(plot_grid, gm$pi, gm$mu, gm$sigma),
                                 trapz_weights(plot_grid))
  refit_components <- sapply(seq_len(refit$Khat), function(k)
    refit$pi[k] * daepd(plot_grid, refit$mu[k], refit$sigma[k],
                        refit$alpha[k], refit$tau[k]))
  colnames(refit_components) <- sprintf("refit_comp_%d", seq_len(refit$Khat))
  kde_plot <- approx(kde$grid, kde$ghat, xout = plot_grid, rule = 2)$y
  ymax <- max(f_raw, f_refit, f_gm, kde_plot, refit_components,
              hist(x, breaks = 25, plot = FALSE)$density) * 1.05

  components_file <- file.path(output_dir, "galaxies_pmhd_refit_components.csv")
  write.csv(data.frame(velocity = plot_grid,
                       refit_total = f_refit,
                       refit_components,
                       check.names = FALSE),
            components_file, row.names = FALSE)

  plot_file <- file.path(output_dir, "galaxies_pmhd_refit_fit.pdf")
  pdf(plot_file, width = 9, height = 5.5)
  par(mar = c(4.2, 4.2, 1.2, 0.8))
  hist(x, breaks = 25, freq = FALSE,
       col = "grey90", border = "grey60",
       xlim = range(plot_grid), ylim = c(0, ymax),
       xlab = "velocity (km/s)", ylab = "density", main = "")
  rug(x, col = "grey20", lwd = 0.6)
  lines(plot_grid, kde_plot, col = "grey40", lty = 3, lwd = 1.2)
  lines(plot_grid, f_raw, col = "red3", lwd = 1.6, lty = 2)
  comp_cols <- grDevices::hcl.colors(refit$Khat, palette = "Dark 3")
  for (k in seq_len(refit$Khat)) {
    lines(plot_grid, refit_components[, k],
          col = comp_cols[k], lty = 3, lwd = 1.05)
  }
  lines(plot_grid, f_refit, col = "firebrick", lwd = 2.2)
  lines(plot_grid, f_gm, col = "navy", lwd = 1.4, lty = 4)
  legend("topright",
         legend = c(
           "histogram", "KDE",
           sprintf("PMHD-MCP (K=%d, H=%.4f)",
                   raw_fit$Khat, raw_H_check),
           sprintf("PMHD-MCP-R (K=%d, H=%.4f)",
                   refit$Khat, refit_H),
           sprintf("Gaussian mixture (K=%d, H=%.4f)",
                   gm$Khat, gm_H),
           "PMHD refit components"
         ),
         pch = c(22, NA, NA, NA, NA, NA),
         pt.bg = c("grey90", NA, NA, NA, NA, NA),
         pt.cex = c(1.6, NA, NA, NA, NA, NA),
         lty = c(NA, 3, 2, 1, 4, 3),
         col = c("grey60", "grey40", "red3", "firebrick", "navy", "grey25"),
         lwd = c(NA, 1.2, 1.6, 2.2, 1.4, 1.05),
         bty = "n", cex = 0.85)
  dev.off()

  out <- list(
    pmhd = pmhd,
    raw_fit = raw_fit,
    refit = refit,
    gm = gm,
    summary = summary_df,
    table = table_out,
    refit_components = data.frame(velocity = plot_grid,
                                  refit_total = f_refit,
                                  refit_components,
                                  check.names = FALSE),
    kde = kde,
    config = list(
      K_max = K_max,
      lambda_grid = lambda_grid,
      h_elbow_jump_frac = h_elbow_jump_frac,
      h_elbow_jump_min = h_elbow_jump_min,
      gamma_anchor = gamma_anchor,
      shrink_floor = shrink_floor,
      raw_trust = raw_trust,
      fast_paper_mode = fast_paper_mode,
      fixed_lambda = fixed_lambda,
      refit_time = refit_time
    )
  )
  rds_file <- file.path(output_dir, "galaxies_pmhd_refit_results.rds")
  saveRDS(out, rds_file)

  cat("\nPMHD refit summary:\n")
  print(summary_df)
  cat("\nFiles written:\n")
  cat("  ", table_file, "\n", sep = "")
  cat("  ", summary_file, "\n", sep = "")
  cat("  ", components_file, "\n", sep = "")
  cat("  ", plot_file, "\n", sep = "")
  cat("  ", rds_file, "\n", sep = "")
  invisible(out)
}

run_galaxies_final_fast <- function(
    output_dir = file.path(galaxies_dir, "K6fast"),
    source_dirs = c(file.path(galaxies_dir, "K6final"),
                    file.path(galaxies_dir, "K6"))) {
  sync_galaxies_final_outputs(output_dir = output_dir,
                              source_dirs = source_dirs)
}

if (identical(environment(), globalenv())) {
  fast_flag <- tolower(Sys.getenv("GALAXIES_FAST_MODE", unset = ""))
  if (fast_flag %in% c("1", "true", "yes", "y")) {
    run_galaxies_final_fast()
  } else {
    run_galaxies_final_analysis()
  }
}

