#!/usr/bin/env Rscript

# Reproducible density analysis of the corrected MASS galaxies data.
#
# The default target is the Gaussian KDE with the Sheather--Jones direct
# plug-in bandwidth used in the Venables--Ripley MASS example.  The option
# --bandwidth=nrd0-0.8 instead reproduces the Epanechnikov KDE and bandwidth
# rule used in Scenarios 1 and 2.
#
# Main paper analysis:
#   Rscript galaxies_analysis.R
#
# Quick smoke test:
#   Rscript galaxies_analysis.R --quick --force
#
# The lambda path is checkpointed after every fitted value and can be resumed
# after interruption. Order selection remains part of the fitting run rather
# than a separate post-processing script.

`%||%` <- function(x, y) if (is.null(x)) y else x

.galaxies_script_file <- local({
  arguments <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", arguments, value = TRUE)
  candidate <- if (length(hit))
    sub("^--file=", "", hit[1L]) else
    tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(candidate) || !nzchar(candidate))
    candidate <- file.path(getwd(), "galaxies_analysis.R")
  normalizePath(candidate, winslash = "/", mustWork = FALSE)
})
.galaxies_script_dir <- dirname(.galaxies_script_file)
.galaxies_base_file <- file.path(.galaxies_script_dir, "galaxies_base.R")
if (!file.exists(.galaxies_base_file))
  stop("galaxies_base.R must be beside galaxies_analysis.R.")
source(.galaxies_base_file, local = FALSE, chdir = TRUE)

for (.variable in c(
  "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
)) {
  if (!nzchar(Sys.getenv(.variable)))
    do.call(Sys.setenv, setNames(list("1"), .variable))
}
rm(.variable)

read_galaxies <- function() {
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("R package 'MASS' is required.")
  holder <- new.env(parent = emptyenv())
  utils::data("galaxies", package = "MASS", envir = holder)
  if (!exists("galaxies", envir = holder, inherits = FALSE))
    stop("MASS::galaxies could not be loaded.")
  x <- as.numeric(holder$galaxies)
  if (length(x) != 82L)
    stop("MASS::galaxies must contain the 82 observations used by Roeder.")
  if (x[78L] == 26690) {
    # MASS documents 26690 as a typographical error; the corrected value is
    # 26960.  The separate 5607 observation from Postman et al. is not added,
    # because it was omitted by Roeder and is not part of this 82-point set.
    x[78L] <- 26960
  }
  if (x[78L] != 26960 || any(x == 5607))
    stop("Unexpected galaxies data version or correction state.")
  x
}

# ---------------------------------------------------------------------------
# Completed penalized path and trimmed dominant-jump elbow
# ---------------------------------------------------------------------------

isotonic_profile <- function(K, H2) {
  curve <- data.frame(K = as.integer(K), H2 = as.numeric(H2))
  curve <- curve[is.finite(curve$K) & is.finite(curve$H2), , drop = FALSE]
  curve <- curve[!duplicated(curve$K), , drop = FALSE]
  curve <- curve[order(curve$K), , drop = FALSE]
  if (nrow(curve))
    curve$H2 <- -stats::isoreg(curve$K, -curve$H2)$yf
  curve
}

select_trimmed_dominant_jump <- function(compressed) {
  if (is.null(compressed) || nrow(compressed) < 3L)
    return(list(success = FALSE, reason = "fewer than three orders"))
  curve <- isotonic_profile(compressed$K, compressed$H2_star)
  curve <- curve[order(curve$K, decreasing = TRUE), , drop = FALSE]
  jumps <- diff(curve$H2)
  eligible <- if (length(jumps) >= 3L)
    2L:length(jumps) else seq_along(jumps)
  eligible <- eligible[
    curve$K[eligible + 1L] != 1L &
      is.finite(jumps[eligible]) & jumps[eligible] > 0
  ]
  if (!length(eligible))
    return(list(success = FALSE, reason = "no eligible positive jump"))
  index <- eligible[which.max(jumps[eligible])]
  list(
    success = TRUE, reason = "success",
    K = as.integer(curve$K[index]),
    curve = curve, jumps = jumps, eligible = eligible,
    dominant = index, dominant_jump = jumps[index]
  )
}

resize_penalized_init <- function(init, target_K) {
  parameters <- init[c("pi", "mu", "sigma", "alpha", "tau")]
  parameters$pi <- parameters$pi / sum(parameters$pi)
  while (length(parameters$pi) > target_K) {
    pairs <- utils::combn(seq_along(parameters$pi), 2L)
    scale <- max(stats::median(parameters$sigma), .Machine$double.eps)
    distance <- abs(
      parameters$mu[pairs[1L, ]] - parameters$mu[pairs[2L, ]]
    ) / scale
    pair <- pairs[, which.min(distance)]
    a <- pair[1L]
    b <- pair[2L]
    weight <- parameters$pi[c(a, b)]
    total <- sum(weight)
    new_mu <- sum(weight * parameters$mu[c(a, b)]) / total
    new_sigma <- sqrt(sum(weight * (
      parameters$sigma[c(a, b)]^2 +
        (parameters$mu[c(a, b)] - new_mu)^2
    )) / total)
    parameters$pi[a] <- total
    parameters$mu[a] <- new_mu
    parameters$sigma[a] <- new_sigma
    parameters$alpha[a] <-
      sum(weight * parameters$alpha[c(a, b)]) / total
    parameters$tau[a] <- sum(weight * parameters$tau[c(a, b)]) / total
    for (name in names(parameters))
      parameters[[name]] <- parameters[[name]][-b]
  }
  while (length(parameters$pi) < target_K) {
    j <- which.max(parameters$pi * parameters$sigma)
    offset <- 0.25 * parameters$sigma[j]
    for (name in names(parameters))
      parameters[[name]] <- append(
        parameters[[name]], parameters[[name]][j], after = j
      )
    parameters$pi[c(j, j + 1L)] <- parameters$pi[j] / 2
    parameters$mu[j] <- parameters$mu[j] - offset
    parameters$mu[j + 1L] <- parameters$mu[j + 1L] + offset
    parameters$sigma[c(j, j + 1L)] <-
      0.85 * parameters$sigma[c(j, j + 1L)]
  }
  parameters$pi <- parameters$pi / sum(parameters$pi)
  sort_parameters(parameters)
}

fit_from_record <- function(record) {
  if (is.null(record) || !isTRUE(record$success) || is.null(record$fit))
    return(NULL)
  active_part(record$fit)
}

support_gap_requests <- function(
    path_result, fractions = c(0.25, 0.50, 0.75)) {
  path <- path_result$path
  successful <- path[
    path$success & is.finite(path$K_monotone), , drop = FALSE
  ]
  successful <- successful[order(successful$lambda), , drop = FALSE]
  if (nrow(successful) < 2L) return(list())
  requests <- list()
  for (j in which(abs(diff(successful$K_monotone)) > 1L)) {
    left <- successful[j, , drop = FALSE]
    right <- successful[j + 1L, , drop = FALSE]
    high <- max(left$K_monotone, right$K_monotone)
    low <- min(left$K_monotone, right$K_monotone)
    missing <- seq.int(low + 1L, high - 1L)
    lambdas <- if (left$lambda > 0 && right$lambda > 0) {
      exp(
        (1 - fractions) * log(left$lambda) +
          fractions * log(right$lambda)
      )
    } else {
      left$lambda + fractions * (right$lambda - left$lambda)
    }
    left_fit <- fit_from_record(
      path_result$records[[as.character(left$key)]]
    )
    right_fit <- fit_from_record(
      path_result$records[[as.character(right$key)]]
    )
    for (K in missing) {
      requests[[length(requests) + 1L]] <- list(
        K = as.integer(K), lambdas = lambdas,
        warm = Filter(Negate(is.null), list(
          if (is.null(left_fit)) NULL else resize_penalized_init(left_fit, K),
          if (is.null(right_fit)) NULL else resize_penalized_init(right_fit, K)
        ))
      )
    }
  }
  requests
}

best_path_fit_for_order <- function(path_result, K) {
  candidates <- list()
  losses <- numeric()
  for (record in path_result$records) {
    if (isTRUE(record$success) && !is.null(record$fit) &&
        isTRUE(as.integer(record$fit$Khat) == as.integer(K)) &&
        is.finite(record$fit$H2)) {
      candidates[[length(candidates) + 1L]] <- record$fit
      losses <- c(losses, record$fit$H2)
    }
  }
  if (!length(candidates)) return(NULL)
  candidates[[which.min(losses)]]
}

polish_augmented_orders <- function(
    x, path_result, augmented, additions, config, seed) {
  orders <- sort(unique(as.integer(augmented$K)))
  source_fits <- Filter(
    Negate(is.null),
    c(
      lapply(orders, function(K) best_path_fit_for_order(path_result, K)),
      additions
    )
  )
  source_fits <- lapply(
    source_fits, active_part, support_tol = config$support_tol
  )
  polished <- list()
  for (j in seq_along(orders)) {
    K <- orders[j]
    row <- augmented[augmented$K == K, , drop = FALSE][1L, ]
    lambda <- as.numeric(row$lambda_star)
    if (!is.finite(lambda) || lambda <= 0) {
      positive <- path_result$path$lambda[
        path_result$path$success & path_result$path$lambda > 0
      ]
      lambda <- if (length(positive)) min(positive) else 1e-4
    }
    warm <- lapply(source_fits, resize_penalized_init, target_K = K)
    progress_message(
      "Order polishing ", j, "/", length(orders),
      ": K=", K, ", lambda=", format(lambda, digits = 8),
      ", warm starts=", length(warm), "."
    )
    lambdas <- unique(pmax(
      min(config$lambda_grid[config$lambda_grid > 0]),
      lambda * config$polish_lambda_factors
    ))
    for (ell in seq_along(lambdas)) {
      trial_lambda <- lambdas[ell]
      trial <- fit_pmhd_at_lambda(
        x = x, K = K, lambda = trial_lambda,
        kde = path_result$kde,
        nstart = config$polish_nstart,
        gamma = config$gamma, warm = warm,
        seed = seed + 700000L + 1009L * K + 37L * ell,
        fit_control = config$fit_control
      )
      if (isTRUE(trial$success) &&
          as.integer(trial$fit$K_exact) == K) {
        polished[[length(polished) + 1L]] <- trial$fit
        hit <- which(augmented$K == K)[1L]
        if (trial$fit$H2 < augmented$H2_star[hit]) {
          augmented$H2_star[hit] <- trial$fit$H2
          augmented$lambda_star[hit] <- trial$fit$lambda
        }
        progress_message(
          "Order polishing accepted: K=", K,
          ", lambda=", format(trial_lambda, digits = 8),
          ", H2=", format(trial$fit$H2, digits = 7)
        )
      } else {
        attained <- if (isTRUE(trial$success))
          trial$fit$K_exact else "failed"
        progress_message(
          "Order polishing candidate rejected: lambda=",
          format(trial_lambda, digits = 8),
          ", attained K=", attained, "."
        )
      }
    }
  }
  list(augmented = augmented, fits = polished)
}

complete_and_select_pmhd <- function(x, path_result, config, seed) {
  requests <- support_gap_requests(path_result, config$gap_fractions)
  additions <- list()
  progress_message(
    "Gap inspection completed: ", length(requests),
    " missing-order request(s)."
  )
  for (r in seq_along(requests)) {
    request <- requests[[r]]
    candidates <- list()
    progress_message(
      "Gap request ", r, "/", length(requests),
      ": target K=", request$K, ", ",
      length(request$lambdas), " penalized candidate(s)."
    )
    for (j in seq_along(request$lambdas)) {
      lambda <- request$lambdas[j]
      progress_message(
        "Gap request ", r, "/", length(requests),
        ", candidate ", j, "/", length(request$lambdas),
        ": lambda=", format(lambda, digits = 8)
      )
      result <- fit_pmhd_at_lambda(
        x = x, K = request$K, lambda = lambda,
        kde = path_result$kde, nstart = config$gap_nstart,
        gamma = config$gamma, warm = request$warm,
        seed = seed + r * 100003L + as.integer(lambda * 1e6),
        fit_control = config$fit_control
      )
      if (isTRUE(result$success) &&
          isTRUE(as.integer(result$fit$K_exact) == request$K)) {
        candidates[[length(candidates) + 1L]] <- result$fit
        progress_message(
          "Gap candidate accepted: K=", result$fit$K_exact,
          ", H2=", format(result$fit$H2, digits = 7)
        )
      } else {
        attained <- if (isTRUE(result$success))
          result$fit$K_exact else "failed"
        progress_message(
          "Gap candidate rejected: attained K=", attained,
          ", required K=", request$K
        )
      }
    }
    if (length(candidates)) {
      losses <- vapply(candidates, function(fit) fit$H2, numeric(1))
      additions[[length(additions) + 1L]] <-
        candidates[[which.min(losses)]]
    }
  }

  augmented <- path_result$compressed
  for (fit in additions) {
    K <- as.integer(fit$Khat)
    hit <- which(augmented$K == K)
    if (length(hit)) {
      if (fit$H2 < augmented$H2_star[hit[1L]]) {
        augmented$H2_star[hit[1L]] <- fit$H2
        augmented$lambda_star[hit[1L]] <- fit$lambda
      }
    } else {
      row <- augmented[1L, , drop = FALSE]
      row[] <- NA
      row$K <- K
      row$H2_star <- fit$H2
      row$lambda_star <- fit$lambda
      augmented <- rbind(augmented, row)
    }
  }
  augmented <- augmented[order(augmented$K, decreasing = TRUE), , drop = FALSE]
  polishing <- polish_augmented_orders(
    x, path_result, augmented, additions, config, seed
  )
  augmented <- polishing$augmented
  selection <- select_trimmed_dominant_jump(augmented)
  if (!isTRUE(selection$success))
    stop("Trimmed dominant-jump elbow failed: ", selection$reason)
  K <- selection$K
  candidates <- Filter(
    Negate(is.null), list(best_path_fit_for_order(path_result, K))
  )
  candidates <- c(
    candidates,
    Filter(function(fit) as.integer(fit$Khat) == K, additions),
    Filter(function(fit) as.integer(fit$Khat) == K, polishing$fits)
  )
  if (!length(candidates))
    stop("Selected order has no penalized fit.")
  losses <- vapply(candidates, function(fit) fit$H2, numeric(1))
  fit <- candidates[[which.min(losses)]]
  progress_message(
    "Order selection completed: K=", K,
    ", lambda=", format(fit$lambda, digits = 8),
    ", H2=", format(fit$H2, digits = 7)
  )
  list(
    fit = fit, K = K, lambda = fit$lambda,
    triggered = length(requests) > 0L,
    additions = additions, polished = polishing$fits,
    augmented = augmented,
    selection = selection
  )
}

build_display_order_profile <- function(x, result, config, seed) {
  orders <- seq.int(2L, config$K_max)
  rows <- vector("list", length(orders))
  fits <- vector("list", length(orders))
  names(fits) <- as.character(orders)

  for (j in seq_along(orders)) {
    K <- orders[j]
    path_fit <- best_path_fit_for_order(result$path, K)
    candidates <- Filter(
      Negate(is.null),
      c(
        list(path_fit),
        Filter(
          function(fit)
            !is.null(fit) && as.integer(fit$Khat) == K,
          result$selected$additions
        ),
        Filter(
          function(fit)
            !is.null(fit) && as.integer(fit$Khat) == K,
          result$selected$polished
        )
      )
    )
    source <- if (length(candidates)) "lambda path" else
      "fixed-order penalized completion"

    if (!length(candidates)) {
      completion_lambdas <- unique(
        result$selected$lambda * c(0.25, 0.50, 1)
      )
      progress_message(
        "Display-profile completion ", j, "/", length(orders),
        ": fitting missing K=", K, " at ",
        length(completion_lambdas), " penalized lambda value(s)."
      )
      for (ell in seq_along(completion_lambdas)) {
        lambda <- completion_lambdas[ell]
        trial <- fit_pmhd_at_lambda(
          x = x, K = K, lambda = lambda,
          kde = result$path$kde,
          nstart = max(10L, config$gap_nstart),
          gamma = config$gamma,
          warm = list(resize_penalized_init(result$pmhd, K)),
          seed = seed + 900000L + K * 101L + ell,
          fit_control = config$fit_control
        )
        if (isTRUE(trial$success) &&
            as.integer(trial$fit$K_exact) == K)
          candidates[[length(candidates) + 1L]] <- trial$fit
      }
    }
    if (!length(candidates))
      stop("Could not construct the display profile at K=", K, ".")

    losses <- vapply(candidates, function(fit) fit$H2, numeric(1))
    fit <- candidates[[which.min(losses)]]
    fits[[as.character(K)]] <- fit
    rows[[j]] <- data.frame(
      K = K, H2 = fit$H2, H = sqrt(pmax(fit$H2, 0)),
      log_H = log(sqrt(pmax(fit$H2, .Machine$double.xmin))),
      lambda = fit$lambda, source = source,
      stringsAsFactors = FALSE
    )
  }
  list(profile = do.call(rbind, rows), fits = fits)
}

# ---------------------------------------------------------------------------
# Fitting, summaries, and figures
# ---------------------------------------------------------------------------

sort_fit <- function(fit) {
  ordering <- order(fit$mu)
  for (name in c("pi", "mu", "sigma", "alpha", "tau"))
    fit[[name]] <- fit[[name]][ordering]
  fit$Khat <- length(fit$pi)
  fit
}

fit_density <- function(grid, fit) {
  mixture_aepd_grid(
    grid, fit$pi, fit$mu, fit$sigma, fit$alpha, fit$tau
  )
}

fit_gaussian_mixture <- function(x, K_max) {
  if (!requireNamespace("mclust", quietly = TRUE))
    stop("R package 'mclust' is required.")
  if (!"package:mclust" %in% search())
    suppressPackageStartupMessages(library("mclust", character.only = TRUE))
  fit <- mclust::Mclust(x, G = seq_len(K_max), verbose = FALSE)
  list(
    Khat = fit$G, model_name = fit$modelName,
    loglik = fit$loglik, fit = fit
  )
}

component_table <- function(fit, n) {
  data.frame(
    component = seq_len(fit$Khat),
    pi = fit$pi, effective_count = n * fit$pi,
    mu = fit$mu, sigma = fit$sigma,
    alpha = fit$alpha, tau = fit$tau
  )
}

analyze_galaxies_once <- function(
    x, config, seed, path_checkpoint_file = NULL,
    resume_path_checkpoint = TRUE) {
  progress_message(
    "Step 1/5: fitting the PMHD-MCP lambda path (seed=", seed, ")."
  )
  path <- fit_pmhd_path(
    x = x, K = config$K_max,
    lambda_grid = config$lambda_grid,
    nstart = config$nstart, gamma = config$gamma,
    grid_n = config$grid_n,
    seed = seed, support_tol = config$support_tol,
    max_grid_expansions = 3L,
    adaptive_nstart = min(5L, config$nstart),
    lambda_upper_cap = config$lambda_max,
    checkpoint_file = path_checkpoint_file,
    checkpoint_id = list(
      version = config$path_version, seed = seed,
      K_max = config$K_max, lambda_grid = config$lambda_grid,
      nstart = config$nstart, gamma = config$gamma,
      grid_n = config$grid_n,
      fit_control = config$fit_control
    ),
    resume_checkpoint = resume_path_checkpoint,
    fit_control = config$fit_control
  )
  if (!isTRUE(path$success))
    stop("PMHD path failed: ", path$reason)
  progress_message(
    "Step 2/5: completing skipped orders and applying the elbow selector."
  )
  selected <- complete_and_select_pmhd(x, path, config, seed)
  progress_message("Step 3/5: extracting the selected penalized fit.")
  pmhd <- sort_fit(active_part(selected$fit, config$support_tol))
  progress_message("Step 4/5: fitting the Gaussian-mixture comparator.")
  gaussian <- fit_gaussian_mixture(x, config$K_max)
  progress_message("Step 5/5: computing fitted densities and distances.")
  kde <- path$kde
  pmhd_density <- normalize_grid_density(
    fit_density(kde$grid, pmhd), kde$w
  )
  gaussian_density <- normalize_grid_density(
    gm_density_grid(kde$grid, gaussian$fit), kde$w
  )
  H2 <- pmax(c(
    hellinger2_grid(kde$phat, pmhd_density, kde$w),
    hellinger2_grid(kde$phat, gaussian_density, kde$w)
  ), 0)
  summary <- data.frame(
    method = c("PMHD-MCP", "Gaussian mixture"),
    Khat = c(pmhd$Khat, gaussian$Khat),
    H2 = H2,
    H = sqrt(H2),
    lambda = c(selected$lambda, NA_real_),
    stringsAsFactors = FALSE
  )
  density_advantage <- data.frame(
    comparison = "PMHD-MCP relative to Gaussian mixture",
    H_reduction = 1 - summary$H[1L] / summary$H[2L],
    H2_reduction = 1 - summary$H2[1L] / summary$H2[2L],
    reference = paste(
      "corrected-data Gaussian KDE with",
      kde$bandwidth_method, "bandwidth"
    )
  )
  list(
    path = path, selected = selected, pmhd = pmhd,
    gaussian = gaussian, summary = summary,
    density_advantage = density_advantage
  )
}

write_selector_figure <- function(result, file, n) {
  path <- result$path$path
  path <- path[path$success & is.finite(path$H2), , drop = FALSE]
  selected <- result$selected
  K_max <- max(result$display_profile$profile$K)
  positive <- path$lambda[path$lambda > 0]
  floor <- if (length(positive)) min(positive) / 2 else 1e-6
  x <- pmax(path$lambda, floor)
  selected_x <- max(selected$lambda, floor)
  grDevices::pdf(file, width = 10, height = 5.4)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 4.7))
  graphics::plot(
    x, n * path$H2, type = "b", log = "x", pch = 16,
    xlab = expression(lambda), ylab = expression(n * H^2),
    col = "black", main = "Hellinger loss path"
  )
  graphics::abline(v = selected_x, lty = 2, col = "firebrick")
  graphics::par(new = TRUE)
  graphics::plot(
    x, path$K_monotone, type = "s", log = "x", axes = FALSE,
    xlab = "", ylab = "", col = "grey45", lwd = 1.8,
    ylim = c(2, K_max)
  )
  graphics::axis(4, at = seq.int(2L, K_max), col.axis = "grey35")
  graphics::mtext("Monotone K", side = 4, line = 2.8, col = "grey35")
  graphics::legend(
    "topleft",
    c(expression(n * H^2), "monotone K", "selected lambda"),
    col = c("black", "grey45", "firebrick"),
    lty = c(1, 1, 2), pch = c(16, NA, NA), bty = "n"
  )

  profile <- result$display_profile$profile
  graphics::plot(
    profile$K, profile$log_H, type = "b", pch = 16,
    xlab = "penalized order K", ylab = expression(log(H[K])),
    xaxt = "n", xlim = c(2, K_max),
    main = "Completed order profile"
  )
  graphics::axis(1, at = seq.int(2L, K_max))
  graphics::abline(v = selected$K, lty = 2, col = "firebrick")
  selected_row <- profile[profile$K == selected$K, , drop = FALSE]
  graphics::points(
    selected_row$K, selected_row$log_H,
    pch = 21, bg = "firebrick", cex = 1.5
  )
}

write_fit_figure <- function(x, result, file) {
  kde <- result$path$kde
  grid <- kde$grid
  pmhd_density <- normalize_grid_density(
    fit_density(grid, result$pmhd), kde$w
  )
  gaussian_density <- normalize_grid_density(
    gm_density_grid(grid, result$gaussian$fit), kde$w
  )
  components <- vapply(seq_len(result$pmhd$Khat), function(k)
    result$pmhd$pi[k] * daepd(
      grid, result$pmhd$mu[k], result$pmhd$sigma[k],
      result$pmhd$alpha[k], result$pmhd$tau[k]
    ), numeric(length(grid)))
  histogram <- graphics::hist(
    x, breaks = 25, plot = FALSE
  )
  y_upper <- 1.10 * max(
    histogram$density, kde$phat, pmhd_density,
    gaussian_density, components, na.rm = TRUE
  )
  grDevices::pdf(file, width = 9, height = 6.2)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(
    histogram, freq = FALSE,
    col = "grey90", border = "grey60", main = "",
    xlim = range(grid), ylim = c(0, y_upper),
    xlab = expression("velocity (km "*s^{-1}*")"),
    ylab = "density"
  )
  graphics::lines(
    grid, kde$phat, col = "grey35", lwd = 1.8, lty = 3
  )
  for (k in seq_len(result$pmhd$Khat))
    graphics::lines(
      grid, components[, k],
      col = grDevices::adjustcolor("firebrick", 0.55), lty = 3
    )
  graphics::lines(grid, pmhd_density, col = "firebrick", lwd = 2.2)
  graphics::lines(
    grid, gaussian_density, col = "navy", lwd = 2, lty = 2
  )
  graphics::rug(x, col = grDevices::adjustcolor("black", 0.35))
  pmhd_H <- result$summary$H[result$summary$method == "PMHD-MCP"]
  gaussian_H <- result$summary$H[
    result$summary$method == "Gaussian mixture"
  ]
  graphics::legend(
    "topright", inset = 0.01,
    c(
      "histogram", "KDE",
      sprintf(
        "PMHD-MCP (K=%d, H=%.4f)", result$pmhd$Khat, pmhd_H
      ),
      sprintf(
        "Gaussian mixture (K=%d, H=%.4f)",
        result$gaussian$Khat, gaussian_H
      ),
      "PMHD components"
    ),
    col = c("grey60", "grey35", "firebrick", "navy", "firebrick"),
    lty = c(1, 3, 1, 2, 3),
    lwd = c(6, 1.8, 2.2, 2, 1), bty = "n"
  )
  curves <- data.frame(
    velocity = grid, kde = kde$phat,
    pmhd_mcp = pmhd_density, gaussian = gaussian_density
  )
  for (k in seq_len(result$pmhd$Khat))
    curves[[paste0("component_", k)]] <- components[, k]
  curves
}

write_pdf_with_fallback <- function(writer, file, ...) {
  tryCatch(
    {
      value <- writer(..., file = file)
      list(value = value, file = file)
    },
    error = function(error) {
      fallback <- sub("\\.pdf$", "_updated.pdf", file)
      progress_message(
        "Could not overwrite ", basename(file), " (",
        conditionMessage(error), "). Writing ", basename(fallback),
        " instead."
      )
      value <- writer(..., file = fallback)
      list(value = value, file = fallback)
    }
  )
}

write_outputs <- function(result, output_dir, metadata, x) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- result$path$path
  augmented <- result$selected$augmented
  curve <- result$selected$selection$curve
  utils::write.csv(
    data.frame(index = seq_along(x), velocity = x),
    file.path(output_dir, "galaxies_data_used.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(
      kernel = result$path$kde$kernel,
      bandwidth_selector = result$path$kde$bandwidth_method,
      bandwidth_km_s = result$path$kde$bandwidth,
      grid_n = result$path$kde$grid_n,
      observation_78 = x[78L],
      includes_5607 = any(x == 5607)
    ),
    file.path(output_dir, "galaxies_kde_spec.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    path, file.path(output_dir, "galaxies_selector_path.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    augmented, file.path(output_dir, "galaxies_completed_path.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    curve, file.path(output_dir, "galaxies_isotonic_path.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    result$display_profile$profile,
    file.path(output_dir, "galaxies_display_order_profile.csv"),
    row.names = FALSE
  )
  selector_summary <- data.frame(
    selector = "trimmed dominant-jump elbow",
    Khat = result$selected$K,
    lambda = result$selected$lambda,
    dominant_jump = result$selected$selection$dominant_jump,
    gap_completion_triggered = result$selected$triggered,
    completed_orders = paste(
      vapply(result$selected$additions, function(fit) fit$Khat, integer(1)),
      collapse = ","
    ),
    polished_orders = paste(
      vapply(result$selected$polished, function(fit) fit$Khat, integer(1)),
      collapse = ","
    ),
    exact_zero_count = metadata$K_max - result$selected$K
  )
  utils::write.csv(
    selector_summary,
    file.path(output_dir, "galaxies_selector_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    component_table(result$pmhd, length(x)),
    file.path(output_dir, "galaxies_pmhd_components.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    result$summary,
    file.path(output_dir, "galaxies_fit_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    result$density_advantage,
    file.path(output_dir, "galaxies_density_advantage.csv"),
    row.names = FALSE
  )
  selector_pdf <- write_pdf_with_fallback(
    write_selector_figure,
    file.path(output_dir, "galaxies_selector_paths.pdf"),
    result = result, n = length(x)
  )
  fit_pdf <- write_pdf_with_fallback(
    write_fit_figure,
    file.path(output_dir, "galaxies_fit.pdf"),
    x = x, result = result
  )
  curves <- fit_pdf$value
  utils::write.csv(
    curves, file.path(output_dir, "galaxies_fitted_curves.csv"),
    row.names = FALSE
  )
  saveRDS(
    list(metadata = metadata, data = x, result = result),
    file.path(output_dir, "galaxies_results.rds"),
    compress = "gzip"
  )
  message("Selected K=", result$selected$K,
          ", lambda=", signif(result$selected$lambda, 7))
  print(result$summary)
  print(result$density_advantage)
  print(component_table(result$pmhd, length(x)))
}

run_galaxies_analysis <- function(
    output_dir = file.path(.galaxies_script_dir, "galaxies_results"),
    force = FALSE, quick = FALSE,
    K_max = 8L, lambda_n = 81L, lambda_min = 1e-4,
    lambda_max = 0.5, nstart = 20L, gap_nstart = 6L,
    polish_nstart = 12L,
    gamma = 3, grid_n = 1024L,
    bandwidth_method = "sj-dpi",
    alpha_lower = 1, alpha_upper = 8.5,
    seed = 20250420L) {
  if (quick) {
    lambda_n <- 7L
    nstart <- 6L
    gap_nstart <- 3L
    polish_nstart <- 3L
    grid_n <- 192L
  }
  progress_message("Loading MASS::galaxies data.")
  x <- read_galaxies()
  progress_message(
    "Loaded ", length(x),
    " observations; corrected observation 78 to 26960 and retained ",
    "Roeder's omission of 5607."
  )
  bandwidth_method <- match.arg(
    bandwidth_method,
    c("sj-dpi", "nrd0-0.8", "gaussian-nrd0-0.8")
  )
  old_bandwidth_option <- getOption("pmhd.bandwidth_method")
  options(pmhd.bandwidth_method = bandwidth_method)
  on.exit(
    options(pmhd.bandwidth_method = old_bandwidth_option),
    add = TRUE
  )
  bandwidth_tag <- switch(
    bandwidth_method,
    "sj-dpi" = "sjdpi",
    "nrd0-0.8" = "epan_nrd08",
    "gaussian-nrd0-0.8" = "gauss_nrd08"
  )
  bandwidth_label <- if (bandwidth_method == "sj-dpi")
    "Sheather-Jones direct plug-in" else "0.8 * bw.nrd0"
  kernel_label <- if (bandwidth_method == "nrd0-0.8")
    "Epanechnikov" else "Gaussian"
  progress_message(
    "Using ", kernel_label, " KDE with bandwidth ", bandwidth_label, "."
  )
  lambda_grid <- c(
    0, exp(seq(
      log(lambda_min), log(lambda_max), length.out = lambda_n - 1L
    ))
  )
  fit_control <- list(
    alpha_lower = alpha_lower, alpha_upper = alpha_upper,
    tau_lower = 0.01, tau_upper = 0.99,
    max_outer = if (quick) 20L else 40L,
    theta_stage1 = if (quick) 15L else 30L,
    theta_stage2 = if (quick) 30L else 60L,
    weight_max_iter = if (quick) 150L else 300L,
    weight_lla_max_iter = if (quick) 5L else 10L,
    parameter_tol = 1e-2,
    objective_tol = if (quick) 1e-5 else 2e-6,
    weight_pg_tol = if (quick) 1e-5 else 1e-6,
    kkt_tol = if (quick) 1e-4 else 1e-5,
    support_tol = 1e-8
  )
  config <- list(
    version = paste0("galaxies_", bandwidth_tag, "_density_dominant_v2"),
    path_version = paste0("galaxies_", bandwidth_tag, "_density_v1"),
    K_max = K_max, lambda_grid = lambda_grid,
    lambda_max = lambda_max,
    nstart = nstart, gap_nstart = gap_nstart,
    polish_nstart = polish_nstart,
    gamma = gamma, grid_n = grid_n,
    bandwidth_method = bandwidth_method,
    bandwidth_label = bandwidth_label,
    kernel_label = kernel_label,
    alpha_lower = alpha_lower, alpha_upper = alpha_upper, seed = seed,
    support_tol = 1e-8,
    gap_fractions = c(0.25, 0.50, 0.75),
    polish_lambda_factors = c(0.50, 1),
    fit_control = fit_control
  )
  metadata <- c(
    config,
    list(
      n = length(x), data = "MASS::galaxies",
      data_correction = "observation 78: 26690 -> 26960",
      omitted_observation = "5607 (following Roeder 1990)",
      kernel = kernel_label,
      bandwidth = bandwidth_label,
      support = "exact simplex zeros; no pruning threshold"
    )
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_file <- file.path(
    output_dir,
    paste0(
      "galaxies_analysis_checkpoint_", bandwidth_tag,
      "_dominant_v2.rds"
    )
  )
  path_checkpoint_file <- file.path(
    output_dir,
    paste0("galaxies_path_checkpoint_", bandwidth_tag, "_v1.rds")
  )
  result <- NULL
  if (!force && file.exists(checkpoint_file)) {
    progress_message("Reading completed checkpoint: ", checkpoint_file)
    checkpoint <- readRDS(checkpoint_file)
    checkpoint_metadata <- checkpoint$metadata
    for (name in c(
      "bandwidth_method", "bandwidth_label", "kernel_label"
    )) {
      if (is.null(checkpoint_metadata[[name]]))
        checkpoint_metadata[[name]] <- metadata[[name]]
    }
    if (!isTRUE(all.equal(
      checkpoint_metadata, metadata, tolerance = 0,
      check.attributes = FALSE
    ))) {
      stop(
        "Checkpoint settings differ. Use --force or another --out directory."
      )
    }
    result <- checkpoint$result
    progress_message("Completed checkpoint loaded.")
  } else {
    progress_message(
      "No reusable completed checkpoint; starting the full analysis."
    )
    result <- analyze_galaxies_once(
      x, config, seed,
      path_checkpoint_file = path_checkpoint_file,
      resume_path_checkpoint = !force
    )
  }
  if (is.null(result$display_profile)) {
    progress_message(
      "Constructing the K=2,...,", K_max,
      " penalized display profile."
    )
    result$display_profile <- build_display_order_profile(
      x, result, config, seed
    )
  }
  progress_message("Saving completed checkpoint: ", checkpoint_file)
  saveRDS(
    list(metadata = metadata, result = result),
    checkpoint_file, compress = "gzip"
  )
  progress_message("Writing tables, figures, and the result object.")
  write_outputs(result, output_dir, metadata, x)
  progress_message("Main galaxies analysis completed.")

  invisible(result)
}

parse_cli <- function(arguments) {
  config <- list(
    output_dir = file.path(.galaxies_script_dir, "galaxies_results"),
    force = FALSE, quick = FALSE,
    K_max = 8L, lambda_n = 81L, lambda_min = 1e-4,
    lambda_max = 0.5, nstart = 20L, gap_nstart = 6L,
    polish_nstart = 12L,
    grid_n = 1024L, bandwidth_method = "sj-dpi",
    seed = 20250420L
  )
  value_after <- function(argument, prefix)
    sub(paste0("^", prefix), "", argument)
  for (argument in arguments) {
    if (argument == "--force") config$force <- TRUE
    else if (argument == "--quick") config$quick <- TRUE
    else if (grepl("^--out=", argument))
      config$output_dir <- value_after(argument, "--out=")
    else if (grepl("^--K-max=", argument))
      config$K_max <- as.integer(value_after(argument, "--K-max="))
    else if (grepl("^--lambda-n=", argument))
      config$lambda_n <- as.integer(value_after(argument, "--lambda-n="))
    else if (grepl("^--lambda-min=", argument))
      config$lambda_min <- as.numeric(value_after(argument, "--lambda-min="))
    else if (grepl("^--lambda-max=", argument))
      config$lambda_max <- as.numeric(value_after(argument, "--lambda-max="))
    else if (grepl("^--nstart=", argument))
      config$nstart <- as.integer(value_after(argument, "--nstart="))
    else if (grepl("^--gap-nstart=", argument))
      config$gap_nstart <- as.integer(
        value_after(argument, "--gap-nstart=")
      )
    else if (grepl("^--polish-nstart=", argument))
      config$polish_nstart <- as.integer(
        value_after(argument, "--polish-nstart=")
      )
    else if (grepl("^--grid=", argument))
      config$grid_n <- as.integer(value_after(argument, "--grid="))
    else if (grepl("^--bandwidth=", argument))
      config$bandwidth_method <- value_after(argument, "--bandwidth=")
    else if (grepl("^--seed=", argument))
      config$seed <- as.integer(value_after(argument, "--seed="))
    else stop("Unknown option: ", argument)
  }
  config
}

if (sys.nframe() == 0L)
  do.call(
    run_galaxies_analysis,
    parse_cli(commandArgs(trailingOnly = TRUE))
  )
