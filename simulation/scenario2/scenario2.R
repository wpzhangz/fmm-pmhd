#!/usr/bin/env Rscript

# Reproducible Scenario 2 simulation for Table 2.
#
# Full run (Linux/Windows):
#   Rscript scenario2.R --B=100 --cores=8
#
# Smoke test:
#   Rscript scenario2.R --quick --force
#
# Every replication performs the complete analysis in one pass:
#   clean Scenario 1(b) generation -> Cauchy gross-error contamination ->
#   PMHD-MCP path -> penalized support-gap completion -> trimmed
#   parsimonious elbow -> metrics, together with GM-MLE and AEPD-MLE-HQ
#   on exactly the same contaminated data.
#
# The cache contains completed replication results only. It is used solely
# for checkpoint/resume and is never an input to a separate tuning step.

`%||%` <- function(x, y) if (is.null(x)) y else x

.scenario2_script_file <- local({
  command <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", command, value = TRUE)
  candidate <- if (length(hit)) sub("^--file=", "", hit[1L]) else {
    tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  }
  if (is.null(candidate) || !nzchar(candidate))
    candidate <- file.path(getwd(), "scenario2.R")
  normalizePath(candidate, winslash = "/", mustWork = FALSE)
})
.scenario2_script_dir <- dirname(.scenario2_script_file)
.scenario2_base_file <- file.path(.scenario2_script_dir, "scenario2_base.R")
if (!file.exists(.scenario2_base_file))
  stop("scenario2_base.R must be beside scenario2.R.")
source(.scenario2_base_file, local = FALSE, chdir = TRUE)

for (.variable in c(
  "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
)) {
  if (!nzchar(Sys.getenv(.variable)))
    do.call(Sys.setenv, setNames(list("1"), .variable))
}
rm(.variable)

scenario2_truth <- list(
  label = "Scenario 2", seed = 20201L,
  pi = c(0.3, 0.4, 0.3), mu = c(-4, 0, 4),
  sigma = c(0.7, 0.9, 0.7),
  alpha = c(1.8, 1.7, 1.8), tau = c(0.45, 0.50, 0.55)
)

scenario2_design <- data.frame(
  setting = c(
    "n100_e005", "n100_e010", "n100_e015",
    "n200_e005", "n200_e010", "n200_e015",
    "n500_e020"
  ),
  n = c(100L, 100L, 100L, 200L, 200L, 200L, 500L),
  epsilon = c(0.05, 0.10, 0.15, 0.05, 0.10, 0.15, 0.20),
  stringsAsFactors = FALSE
)

table2_methods <- c("PMHD-MCP", "GM-MLE", "AEPD-MLE-HQ")

empty_metrics <- function(elapsed = NA_real_, error = "",
                          correct = NA_integer_) {
  list(
    fit_ok = 0L, K = NA_integer_, correct = correct,
    H = NA_real_, ISE = NA_real_,
    mse_mu = NA_real_, mse_sigma = NA_real_,
    mse_alpha = NA_real_, mse_tau = NA_real_,
    elapsed = elapsed, error = error,
    selected_lambda = NA_real_, selector_unresolved = NA_integer_,
    completion_triggered = NA_integer_, completed_orders = NA_integer_,
    diagnostics = NULL
  )
}

metrics_from_aepd <- function(fit, truth, evaluation, elapsed) {
  fitted_density <- mixture_aepd_grid(
    evaluation$grid, fit$pi, fit$mu, fit$sigma, fit$alpha, fit$tau
  )
  fitted_density <- normalize_grid_density(fitted_density, evaluation$w)
  h2 <- hellinger2_grid(
    evaluation$density, fitted_density, evaluation$w
  )
  mse <- parameter_mse(fit, truth)
  list(
    fit_ok = 1L, K = fit$Khat,
    correct = as.integer(fit$Khat == length(truth$pi)),
    H = sqrt(pmax(h2, 0)),
    ISE = ise_grid(evaluation$density, fitted_density, evaluation$w),
    mse_mu = unname(mse["mu"]),
    mse_sigma = unname(mse["sigma"]),
    mse_alpha = unname(mse["alpha"]),
    mse_tau = unname(mse["tau"]),
    elapsed = elapsed, error = "",
    selected_lambda = NA_real_, selector_unresolved = 0L,
    completion_triggered = NA_integer_, completed_orders = NA_integer_,
    diagnostics = NULL
  )
}

# -------------------------------------------------------------------------
# PMHD-MCP selector shared with the updated Table 1 program
# -------------------------------------------------------------------------

isotonic_profile <- function(K, H2) {
  curve <- data.frame(K = as.integer(K), H2 = as.numeric(H2))
  curve <- curve[is.finite(curve$K) & is.finite(curve$H2), , drop = FALSE]
  curve <- curve[!duplicated(curve$K), , drop = FALSE]
  curve <- curve[order(curve$K), , drop = FALSE]
  if (nrow(curve))
    curve$H2 <- -stats::isoreg(curve$K, -curve$H2)$yf
  curve
}

select_trimmed_parsimonious <- function(compressed, eta = 0.35) {
  if (is.null(compressed) || nrow(compressed) < 3L)
    return(list(success = FALSE, reason = "fewer than three orders"))
  curve <- isotonic_profile(compressed$K, compressed$H2_star)
  curve <- curve[order(curve$K, decreasing = TRUE), , drop = FALSE]
  jumps <- diff(curve$H2)

  # Remove the first high-order transition and any transition ending at K=1.
  eligible <- if (length(jumps) >= 3L)
    2L:length(jumps) else seq_along(jumps)
  eligible <- eligible[
    curve$K[eligible + 1L] != 1L &
      is.finite(jumps[eligible]) & jumps[eligible] > 0
  ]
  if (!length(eligible))
    return(list(success = FALSE, reason = "no eligible positive jump"))

  substantial <- eligible[
    jumps[eligible] >= eta * max(jumps[eligible])
  ]
  if (!length(substantial))
    return(list(success = FALSE, reason = "no substantial jump"))

  index <- max(substantial)
  list(
    success = TRUE, reason = "success",
    K = as.integer(curve$K[index]),
    curve = curve, jumps = jumps, eligible = eligible,
    substantial = substantial, eta = eta
  )
}

resize_penalized_init <- function(init, target_K) {
  pars <- init[c("pi", "mu", "sigma", "alpha", "tau")]
  pars$pi <- pars$pi / sum(pars$pi)
  while (length(pars$pi) > target_K) {
    pairs <- utils::combn(seq_along(pars$pi), 2L)
    scale <- max(stats::median(pars$sigma), .Machine$double.eps)
    distance <- abs(pars$mu[pairs[1L, ]] - pars$mu[pairs[2L, ]]) / scale
    pair <- pairs[, which.min(distance)]
    a <- pair[1L]
    b <- pair[2L]
    weight <- pars$pi[c(a, b)]
    total <- sum(weight)
    new_mu <- sum(weight * pars$mu[c(a, b)]) / total
    new_sigma <- sqrt(sum(weight * (
      pars$sigma[c(a, b)]^2 + (pars$mu[c(a, b)] - new_mu)^2
    )) / total)
    pars$pi[a] <- total
    pars$mu[a] <- new_mu
    pars$sigma[a] <- new_sigma
    pars$alpha[a] <- sum(weight * pars$alpha[c(a, b)]) / total
    pars$tau[a] <- sum(weight * pars$tau[c(a, b)]) / total
    for (name in names(pars)) pars[[name]] <- pars[[name]][-b]
  }
  while (length(pars$pi) < target_K) {
    j <- which.max(pars$pi * pars$sigma)
    offset <- 0.25 * pars$sigma[j]
    for (name in names(pars))
      pars[[name]] <- append(pars[[name]], pars[[name]][j], after = j)
    pars$pi[c(j, j + 1L)] <- pars$pi[j] / 2
    pars$mu[j] <- pars$mu[j] - offset
    pars$mu[j + 1L] <- pars$mu[j + 1L] + offset
    pars$sigma[c(j, j + 1L)] <- 0.85 * pars$sigma[c(j, j + 1L)]
  }
  pars$pi <- pars$pi / sum(pars$pi)
  sort_parameters(pars)
}

fit_from_record <- function(record) {
  if (is.null(record) || !isTRUE(record$success) || is.null(record$fit))
    return(NULL)
  active_part(record$fit)
}

support_gap_requests <- function(path_result,
                                 fractions = c(0.25, 0.50, 0.75)) {
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
      exp((1 - fractions) * log(left$lambda) +
            fractions * log(right$lambda))
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
          if (is.null(left_fit)) NULL else
            resize_penalized_init(left_fit, K),
          if (is.null(right_fit)) NULL else
            resize_penalized_init(right_fit, K)
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

complete_and_select_pmhd <- function(x, path_result, config, seed) {
  requests <- support_gap_requests(path_result, config$gap_fractions)
  additions <- list()

  for (r in seq_along(requests)) {
    request <- requests[[r]]
    candidates <- list()
    for (lambda in request$lambdas) {
      result <- fit_pmhd_at_lambda(
        x = x, K = request$K, lambda = lambda,
        kde = path_result$kde, nstart = config$gap_nstart,
        gamma = config$gamma, warm = request$warm,
        seed = seed + r * 100003L + as.integer(lambda * 1e6),
        fit_control = config$pmhd_fit_control
      )
      if (isTRUE(result$success) &&
          isTRUE(as.integer(result$fit$K_exact) == request$K)) {
        candidates[[length(candidates) + 1L]] <- result$fit
      }
    }
    if (length(candidates)) {
      loss <- vapply(candidates, function(fit) fit$H2, numeric(1))
      additions[[length(additions) + 1L]] <-
        candidates[[which.min(loss)]]
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
  selection <- select_trimmed_parsimonious(augmented, config$eta)
  if (!isTRUE(selection$success)) {
    return(list(
      success = FALSE, reason = selection$reason,
      fit = NULL, K = NA_integer_, lambda = NA_real_,
      triggered = length(requests) > 0L,
      additions = additions, augmented = augmented,
      selection = selection
    ))
  }

  K <- selection$K
  candidates <- list(best_path_fit_for_order(path_result, K))
  candidates <- c(
    Filter(Negate(is.null), candidates),
    Filter(function(fit) as.integer(fit$Khat) == K, additions)
  )
  if (!length(candidates)) {
    return(list(
      success = FALSE, reason = "selected order has no penalized fit",
      fit = NULL, K = K, lambda = NA_real_,
      triggered = length(requests) > 0L,
      additions = additions, augmented = augmented,
      selection = selection
    ))
  }
  loss <- vapply(candidates, function(fit) fit$H2, numeric(1))
  fit <- candidates[[which.min(loss)]]
  list(
    success = TRUE, reason = "success",
    fit = fit, K = K, lambda = fit$lambda,
    triggered = length(requests) > 0L,
    additions = additions, augmented = augmented,
    selection = selection
  )
}

fit_pmhd_mcp <- function(x, truth, evaluation, config, seed) {
  started <- proc.time()[["elapsed"]]
  path_error <- ""
  path_result <- tryCatch(fit_pmhd_path(
    x = x, K = config$K_max,
    lambda_grid = config$lambda_grid,
    nstart = config$pmhd_nstart,
    gamma = config$gamma,
    bw_adjust = config$bw_adjust,
    grid_n = config$pmhd_grid_n,
    seed = seed + 1000000L,
    support_tol = config$support_tol,
    max_grid_expansions = config$max_grid_expansions,
    adaptive_nstart = config$adaptive_nstart,
    fit_control = config$pmhd_fit_control
  ), error = function(e) {
    path_error <<- conditionMessage(e)
    NULL
  })
  if (is.null(path_result) || !isTRUE(path_result$success)) {
    reason <- if (nzchar(path_error)) path_error else
      path_result$reason %||% "path generation failed"
    failure <- empty_metrics(
      proc.time()[["elapsed"]] - started, reason, correct = 0L
    )
    failure$selector_unresolved <- 1L
    return(failure)
  }

  selected <- complete_and_select_pmhd(
    x, path_result, config, seed = seed
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if (!isTRUE(selected$success)) {
    failure <- empty_metrics(elapsed, selected$reason, correct = 0L)
    failure$selector_unresolved <- 1L
    failure$completion_triggered <- as.integer(selected$triggered)
    failure$completed_orders <- length(selected$additions)
    return(failure)
  }

  active <- active_part(selected$fit, config$support_tol)
  active$Khat <- selected$K
  metrics <- metrics_from_aepd(active, truth, evaluation, elapsed)
  metrics$selected_lambda <- selected$lambda
  metrics$selector_unresolved <- 0L
  metrics$completion_triggered <- as.integer(selected$triggered)
  metrics$completed_orders <- length(selected$additions)
  metrics$diagnostics <- list(
    selected_K = selected$K,
    selected_lambda = selected$lambda,
    completed_K = vapply(
      selected$additions, function(fit) as.integer(fit$Khat), integer(1)
    )
  )
  metrics
}

# -------------------------------------------------------------------------
# Competitors evaluated on exactly the same generated sample
# -------------------------------------------------------------------------

fit_gaussian_mle <- function(x, truth, evaluation, K_max = 5L) {
  started <- proc.time()[["elapsed"]]
  error_text <- ""
  fit <- tryCatch({
    if (!requireNamespace("mclust", quietly = TRUE))
      stop("Package 'mclust' is required.")
    if (!"package:mclust" %in% search())
      suppressPackageStartupMessages(library("mclust", character.only = TRUE))
    mclust::Mclust(x, G = seq_len(K_max), verbose = FALSE)
  }, error = function(e) {
    error_text <<- conditionMessage(e)
    NULL
  })
  elapsed <- proc.time()[["elapsed"]] - started
  if (is.null(fit)) return(empty_metrics(elapsed, error_text))
  density <- normalize_grid_density(
    gm_density_grid(evaluation$grid, fit), evaluation$w
  )
  mse <- gm_parameter_mse(fit, truth)
  h2 <- hellinger2_grid(evaluation$density, density, evaluation$w)
  list(
    fit_ok = 1L, K = fit$G,
    correct = as.integer(fit$G == length(truth$pi)),
    H = sqrt(pmax(h2, 0)),
    ISE = ise_grid(evaluation$density, density, evaluation$w),
    mse_mu = unname(mse["mu"]), mse_sigma = unname(mse["sigma"]),
    mse_alpha = NA_real_, mse_tau = NA_real_,
    elapsed = elapsed, error = "",
    selected_lambda = NA_real_, selector_unresolved = NA_integer_,
    completion_triggered = NA_integer_, completed_orders = NA_integer_,
    diagnostics = NULL
  )
}

fit_aepd_likelihood <- function(x, truth, evaluation, config, seed) {
  started <- proc.time()[["elapsed"]]
  error_text <- ""
  fit <- tryCatch(fit_aepd_mle_hq(
    x, K_max = config$K_max, nstart = config$aepd_nstart,
    alpha_lower = 1, alpha_upper = 3, seed = seed + 2000000L,
    max_iter = config$aepd_max_iter,
    maxeval = config$aepd_maxeval
  ), error = function(e) {
    error_text <<- conditionMessage(e)
    NULL
  })
  elapsed <- proc.time()[["elapsed"]] - started
  if (is.null(fit))
    return(empty_metrics(
      elapsed,
      if (nzchar(error_text)) error_text else "AEPD-MLE-HQ failed"
    ))
  metrics <- metrics_from_aepd(fit, truth, evaluation, elapsed)
  metrics$diagnostics <- NULL
  metrics
}

one_replication <- function(b, truth, evaluation, config,
                            clean_seed_base, contam_seed_base) {
  clean_seed <- clean_seed_base + as.integer(b)
  set.seed(clean_seed)
  clean <- rmix_aepd(
    config$n, truth$pi, truth$mu, truth$sigma, truth$alpha, truth$tau
  )
  set.seed(contam_seed_base + as.integer(b))
  number_contaminated <- stats::rbinom(1L, config$n, config$epsilon)
  if (number_contaminated > 0L) {
    contamination <- stats::rcauchy(number_contaminated)
    x <- c(
      clean[seq_len(config$n - number_contaminated)],
      contamination
    )
  } else {
    x <- clean
  }
  fit_seed <- contam_seed_base + 1000000L + as.integer(b)
  list(
    "PMHD-MCP" = fit_pmhd_mcp(
      x, truth, evaluation, config, seed = fit_seed
    ),
    "GM-MLE" = fit_gaussian_mle(
      x, truth, evaluation, K_max = config$K_max
    ),
    "AEPD-MLE-HQ" = fit_aepd_likelihood(
      x, truth, evaluation, config, seed = fit_seed
    )
  )
}

# -------------------------------------------------------------------------
# Summaries, Table 2, checkpointing, and CLI
# -------------------------------------------------------------------------

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

summarize_method <- function(replicates, method, setting, n, epsilon) {
  values <- lapply(replicates, `[[`, method)
  number <- function(name) vapply(values, function(value) {
    item <- value[[name]]
    if (is.null(item) || !length(item)) NA_real_
    else suppressWarnings(as.numeric(item[1L]))
  }, numeric(1))
  character_value <- function(name) vapply(values, function(value) {
    item <- value[[name]]
    if (is.null(item) || !length(item) || is.na(item[1L])) ""
    else as.character(item[1L])
  }, character(1))
  error <- character_value("error")
  error <- error[nzchar(error)]
  correct <- number("correct")
  # Any fitting or selection failure counts as incorrect.
  correct[!is.finite(correct)] <- 0
  data.frame(
    setting = setting, scenario = "Scenario 2", n = n,
    epsilon = epsilon, method = method,
    P_K0 = mean(correct == 1, na.rm = TRUE),
    mean_K = finite_mean(number("K")),
    H = finite_mean(number("H")),
    ISE = finite_mean(number("ISE")),
    MSE_mu = finite_mean(number("mse_mu")),
    MSE_sigma = finite_mean(number("mse_sigma")),
    MSE_alpha = finite_mean(number("mse_alpha")),
    MSE_tau = finite_mean(number("mse_tau")),
    N_resolved = sum(is.finite(number("K"))),
    N_MSE = sum(is.finite(number("mse_mu"))),
    failure_rate = mean(number("fit_ok") != 1),
    unresolved_rate = if (method == "PMHD-MCP")
      mean(number("selector_unresolved") == 1, na.rm = TRUE) else NA_real_,
    error_count = length(error),
    first_error = if (length(error)) error[1L] else "",
    stringsAsFactors = FALSE
  )
}

write_table2 <- function(summary, path, B) {
  fmt <- function(x, digits) {
    if (!is.finite(x)) "--" else formatC(x, format = "f", digits = digits)
  }
  lines <- c(
    sprintf("%% Scenario 2 Table 2; B=%d; generated by scenario2.R", B),
    "\\begin{tabular}{lllcccccccc}",
    "\\toprule",
    paste(
      "$n$ & $\\epsilon$ & Method & $P(\\hat K=K_0)$",
      "& $\\bar K$ & $H$ & ISE",
      "& $\\mathrm{MSE}_\\mu$ & $\\mathrm{MSE}_\\sigma$",
      "& $\\mathrm{MSE}_\\alpha$ & $\\mathrm{MSE}_\\tau$ \\\\"
    ),
    "\\midrule"
  )
  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    lines <- c(lines, sprintf(
      "%d & %.3f & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
      row$n, row$epsilon, row$method,
      fmt(row$P_K0, 3), fmt(row$mean_K, 3),
      fmt(row$H, 4), fmt(row$ISE, 4),
      fmt(row$MSE_mu, 4), fmt(row$MSE_sigma, 4),
      fmt(row$MSE_alpha, 4), fmt(row$MSE_tau, 4)
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, path, useBytes = TRUE)
}

run_scenario2 <- function(
    B = 100L, n_values = c(100L, 200L, 500L),
    eps_values = c(0.05, 0.10, 0.15, 0.20), cores = 4L,
    settings = scenario2_design$setting,
    out_dir = file.path(.scenario2_script_dir, "scenario2_results"),
    force = FALSE, checkpoint_every = NULL,
    pmhd_grid_n = 512L, evaluation_grid_n = 4096L,
    pmhd_nstart = 20L, gap_nstart = 6L, aepd_nstart = 20L,
    quick = FALSE) {
  design <- scenario2_design[
    scenario2_design$setting %in% settings &
      scenario2_design$n %in% as.integer(n_values) &
      scenario2_design$epsilon %in% as.numeric(eps_values),
    , drop = FALSE
  ]
  if (!nrow(design)) stop("No Table 2 design cell was requested.")
  if (B < 1L || cores < 1L || any(n_values < 2L))
    stop("B and cores must be positive, and n must be at least two.")
  if (is.null(checkpoint_every))
    checkpoint_every <- max(1L, 2L * min(as.integer(cores), as.integer(B)))

  if (quick) {
    B <- 1L
    design <- scenario2_design[1L, , drop = FALSE]
    cores <- 1L
    pmhd_grid_n <- 192L
    evaluation_grid_n <- 1024L
    pmhd_nstart <- 3L
    gap_nstart <- 3L
    aepd_nstart <- 1L
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(out_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  pmhd_fit_control <- list(
    alpha_lower = 1, alpha_upper = 3,
    tau_lower = 0.01, tau_upper = 0.99,
    max_outer = if (quick) 25L else 40L,
    theta_stage1 = if (quick) 20L else 30L,
    theta_stage2 = if (quick) 40L else 60L,
    weight_max_iter = if (quick) 200L else 300L,
    weight_lla_max_iter = if (quick) 5L else 10L,
    parameter_tol = 1e-2,
    objective_tol = if (quick) 5e-6 else 2e-6,
    weight_pg_tol = if (quick) 1e-5 else 1e-6,
    kkt_tol = if (quick) 1e-4 else 1e-5,
    support_tol = 1e-8
  )
  version <- "scenario2_table2_direct_v1"
  rows <- list()
  truth <- scenario2_truth
  evaluation_grid <- seq(-22, 22, length.out = evaluation_grid_n)
  evaluation_weights <- trapz_weights(evaluation_grid)
  evaluation <- list(
    grid = evaluation_grid, w = evaluation_weights,
    density = normalize_grid_density(
      mixture_aepd_grid(
        evaluation_grid, truth$pi, truth$mu, truth$sigma,
        truth$alpha, truth$tau
      ),
      evaluation_weights
    )
  )

  for (cell_index in seq_len(nrow(design))) {
      setting_name <- design$setting[cell_index]
      n <- as.integer(design$n[cell_index])
      epsilon <- as.numeric(design$epsilon[cell_index])
      n_slot <- match(n, c(100L, 200L, 500L))
      epsilon_slot <- match(epsilon, c(0.05, 0.10, 0.15, 0.20))
      clean_seed_base <- truth$seed + 100000L * n_slot
      contam_seed_base <- clean_seed_base + 10000000L * epsilon_slot
      config <- list(
        n = n, epsilon = epsilon, contamination = "standard Cauchy",
        K_max = 5L,
        lambda_grid = default_lambda_grid(),
        pmhd_nstart = as.integer(pmhd_nstart),
        gap_nstart = as.integer(gap_nstart),
        aepd_nstart = as.integer(aepd_nstart),
        gamma = 3, eta = 0.35,
        gap_fractions = c(0.25, 0.50, 0.75),
        bw_adjust = 0.80,
        pmhd_grid_n = as.integer(pmhd_grid_n),
        support_tol = 1e-8,
        max_grid_expansions = 3L,
        adaptive_nstart = min(5L, as.integer(pmhd_nstart)),
        pmhd_fit_control = pmhd_fit_control,
        aepd_max_iter = 300L, aepd_maxeval = 250L
      )
      metadata <- list(
        version = version, setting = setting_name, n = n,
        epsilon = epsilon, B = B,
        clean_seed_base = clean_seed_base,
        contam_seed_base = contam_seed_base,
        methods = table2_methods,
        truth = truth, config = config,
        evaluation_grid_n = evaluation_grid_n
      )
      cache_file <- file.path(
        cache_dir,
        sprintf(
          "%s_n%d_eps%03d_B%d.rds",
          version, n, as.integer(round(1000 * epsilon)), B
        )
      )
      replicates <- vector("list", B)
      if (!force && file.exists(cache_file)) {
        cache <- readRDS(cache_file)
        if (!isTRUE(all.equal(
          cache$metadata, metadata, tolerance = 0,
          check.attributes = FALSE
        ))) {
          stop(
            "Cache metadata differ: ", cache_file,
            ". Use --force or a different --out directory."
          )
        }
        replicates <- cache$replicates
        length(replicates) <- B
      }

      missing <- which(vapply(replicates, is.null, logical(1)))
      message(sprintf(
        "Scenario 2, n=%d, epsilon=%.3f: resumed %d/%d",
        n, epsilon, B - length(missing), B
      ))
      if (length(missing)) {
        batches <- split(
          missing,
          ceiling(seq_along(missing) / as.integer(checkpoint_every))
        )
        cluster <- NULL
        if (cores > 1L && length(missing) > 1L) {
          cluster <- parallel::makeCluster(min(cores, length(missing)))
          scenario_file <- .scenario2_script_file
          parallel::clusterExport(
            cluster, "scenario_file", envir = environment()
          )
          parallel::clusterEvalQ(
            cluster,
            suppressWarnings(suppressMessages(
              source(scenario_file, local = FALSE, chdir = TRUE)
            ))
          )
        }
        tryCatch({
          for (batch in batches) {
            results <- if (!is.null(cluster) && length(batch) > 1L) {
              parallel::parLapplyLB(
                cluster, batch, one_replication,
                truth = truth, evaluation = evaluation,
                config = config,
                clean_seed_base = clean_seed_base,
                contam_seed_base = contam_seed_base
              )
            } else {
              lapply(
                batch, one_replication,
                truth = truth, evaluation = evaluation,
                config = config,
                clean_seed_base = clean_seed_base,
                contam_seed_base = contam_seed_base
              )
            }
            replicates[batch] <- results
            saveRDS(
              list(metadata = metadata, replicates = replicates),
              cache_file, compress = "gzip"
            )
            message(sprintf(
              "  checkpoint %d/%d",
              sum(!vapply(replicates, is.null, logical(1))), B
            ))
          }
        }, finally = {
          if (!is.null(cluster))
            try(parallel::stopCluster(cluster), silent = TRUE)
        })
      }

      for (method in table2_methods) {
        row <- summarize_method(
          replicates, method, setting_name, n, epsilon
        )
        rows[[length(rows) + 1L]] <- row
        message(sprintf(
          "  %-12s P(K0)=%.3f H=%.4f ISE=%.4f",
          method, row$P_K0, row$H, row$ISE
        ))
      }
      partial <- do.call(rbind, rows)
      utils::write.csv(
        partial, file.path(out_dir, "scenario2_summary_partial.csv"),
        row.names = FALSE, na = ""
      )
  }

  summary <- do.call(rbind, rows)
  summary_file <- file.path(out_dir, "scenario2_summary.csv")
  table_file <- file.path(out_dir, "scenario2_table2.tex")
  utils::write.csv(summary, summary_file, row.names = FALSE, na = "")
  write_table2(summary, table_file, B)
  saveRDS(
    list(
      version = version, truth = scenario2_truth, design = design,
      summary = summary
    ),
    file.path(out_dir, "scenario2_run_manifest.rds")
  )
  message("Written:")
  message("  ", normalizePath(summary_file, winslash = "/"))
  message("  ", normalizePath(table_file, winslash = "/"))
  invisible(summary)
}

parse_cli <- function(arguments) {
  config <- list(
    B = 100L, n_values = c(100L, 200L, 500L), cores = 4L,
    eps_values = c(0.05, 0.10, 0.15, 0.20),
    settings = scenario2_design$setting,
    out_dir = file.path(.scenario2_script_dir, "scenario2_results"),
    force = FALSE, checkpoint_every = NULL,
    pmhd_grid_n = 512L, evaluation_grid_n = 4096L,
    pmhd_nstart = 20L, gap_nstart = 6L, aepd_nstart = 20L,
    quick = FALSE
  )
  value_after <- function(argument, prefix)
    sub(paste0("^", prefix), "", argument)
  for (argument in arguments) {
    if (argument == "--quick") config$quick <- TRUE
    else if (argument == "--force") config$force <- TRUE
    else if (grepl("^--B=", argument))
      config$B <- as.integer(value_after(argument, "--B="))
    else if (grepl("^--n=", argument))
      config$n_values <- as.integer(strsplit(
        value_after(argument, "--n="), ",", fixed = TRUE
      )[[1L]])
    else if (grepl("^--eps=", argument))
      config$eps_values <- as.numeric(strsplit(
        value_after(argument, "--eps="), ",", fixed = TRUE
      )[[1L]])
    else if (grepl("^--cores=", argument))
      config$cores <- as.integer(value_after(argument, "--cores="))
    else if (grepl("^--settings=", argument))
      config$settings <- strsplit(
        value_after(argument, "--settings="), ",", fixed = TRUE
      )[[1L]]
    else if (grepl("^--out=", argument))
      config$out_dir <- value_after(argument, "--out=")
    else if (grepl("^--checkpoint-every=", argument))
      config$checkpoint_every <- as.integer(
        value_after(argument, "--checkpoint-every=")
      )
    else stop("Unknown option: ", argument)
  }
  config
}

if (sys.nframe() == 0L) {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  do.call(run_scenario2, cli)
}
