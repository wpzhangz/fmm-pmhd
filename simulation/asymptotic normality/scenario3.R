#!/usr/bin/env Rscript

# Reproducible Scenario 3 simulation for Table 3.
#
# Full run:
#   Rscript scenario3.R --B=200 --n=500,1000,2000 --cores=8
#
# Smoke test:
#   Rscript scenario3.R --quick --force
#
# The active order K0=3 is fixed, as required for the active-submodel
# Gaussian approximation.  Each replication directly generates the sample,
# computes the fixed-order minimum-Hellinger estimator at h_n=2*n^(-1/3),
# and stores its root-n error.  Checkpoints contain completed replications
# only and are never used for post-hoc estimator selection.

.scenario3_script_file <- local({
  arguments <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", arguments, value = TRUE)
  candidate <- if (length(hit))
    sub("^--file=", "", hit[1L]) else
    tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(candidate) || !nzchar(candidate))
    candidate <- file.path(getwd(), "scenario3.R")
  normalizePath(candidate, winslash = "/", mustWork = FALSE)
})
.scenario3_script_dir <- dirname(.scenario3_script_file)
.scenario3_base_file <- file.path(.scenario3_script_dir, "scenario3_base.R")
if (!file.exists(.scenario3_base_file))
  stop("scenario3_base.R must be beside scenario3.R.")
source(.scenario3_base_file, local = FALSE, chdir = TRUE)

for (.variable in c(
  "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
)) {
  if (!nzchar(Sys.getenv(.variable)))
    do.call(Sys.setenv, setNames(list("1"), .variable))
}
rm(.variable)

scenario3_truth <- list(
  pi = c(0.3, 0.4, 0.3),
  mu = c(-4, 0, 4),
  sigma = c(0.5, 0.5, 0.5),
  alpha = c(1.4, 1.5, 1.3),
  tau = c(0.45, 0.50, 0.55)
)

valid_estimate <- function(estimate) {
  all(is.finite(unlist(estimate))) &&
    all(estimate$pi > 0.02) &&
    all(estimate$sigma > 0.2 & estimate$sigma < 4) &&
    all(estimate$alpha > 0.5 & estimate$alpha < 4) &&
    all(estimate$tau > 0.02 & estimate$tau < 0.98)
}

fit_active_pmhd <- function(
    x, bandwidth, K, grid_n_max, points_per_bandwidth,
    nstart = 10L, tolerance = 1e-5, max_iter = 100L,
    maxeval = 250L, active_weight_min = 0.01) {
  kde_error <- ""
  kde <- tryCatch(
    make_theory_kde(
      x, bandwidth, grid_n_max = grid_n_max,
      points_per_bandwidth = points_per_bandwidth
    ),
    error = function(e) {
      kde_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(kde))
    return(list(success = FALSE, reason = paste("KDE:", kde_error)))

  alpha_grid <- c(2.0, 1.5, 1.2, 1.0)
  best <- NULL
  best_loss <- Inf
  best_active <- NULL
  best_active_loss <- Inf
  errors <- character()
  for (start in seq_len(nstart)) {
    alpha_start <- alpha_grid[
      ((start - 1L) %% length(alpha_grid)) + 1L
    ]
    initial <- if (start == 1L) NULL else tryCatch(
      init_params_kmeans(x, K, alpha_init = alpha_start),
      error = function(e) NULL
    )
    fit_error <- ""
    fit <- tryCatch(
      pmhd_fixed_order_fit_one(
        x, K, kde$grid, kde$w, kde$phat, init = initial,
        tolerance = tolerance, max_iter = max_iter,
        maxeval = maxeval
      ),
      error = function(e) {
        fit_error <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(fit)) {
      if (nzchar(fit_error)) errors <- c(errors, fit_error)
      next
    }
    if (!all(is.finite(c(
      fit$pi, fit$mu, fit$sigma, fit$alpha, fit$tau, fit$H2
    )))) next
    if (fit$H2 < best_loss) {
      best <- fit
      best_loss <- fit$H2
    }
    if (min(fit$pi) >= active_weight_min && fit$H2 < best_active_loss) {
      best_active <- fit
      best_active_loss <- fit$H2
    }
  }
  selected <- if (!is.null(best_active)) best_active else best
  if (is.null(selected)) {
    reason <- if (length(errors)) errors[1L] else "no valid start"
    return(list(success = FALSE, reason = reason))
  }
  estimate <- align_by_location(selected)
  if (!valid_estimate(estimate))
    return(list(success = FALSE, reason = "invalid estimate"))
  list(
    success = TRUE, reason = "", estimate = estimate,
    H2 = selected$H2, grid_n = kde$grid_n,
    grid_step_over_h = kde$resolution
  )
}

one_replication <- function(
    b, n, truth, theta_true, config, seed_base) {
  seed <- seed_base + as.integer(b)
  set.seed(seed)
  x <- rmix_aepd(
    n, truth$pi, truth$mu, truth$sigma, truth$alpha, truth$tau
  )
  fit <- fit_active_pmhd(
    x = x, bandwidth = config$bandwidth_constant * n^(-config$beta),
    K = length(truth$pi),
    grid_n_max = config$kde_grid_n_max,
    points_per_bandwidth = config$kde_points_per_bandwidth,
    nstart = config$nstart, tolerance = config$tolerance,
    max_iter = config$max_iter, maxeval = config$maxeval,
    active_weight_min = config$active_weight_min
  )
  if (!isTRUE(fit$success))
    return(list(ok = FALSE, reason = fit$reason, seed = seed))
  # The optimizer is data driven.  This branch check is only a diagnostic
  # defining the local active-submodel solution summarized by the theorem.
  if (max(abs(fit$estimate$mu - truth$mu)) > 2)
    return(list(ok = FALSE, reason = "wrong local branch", seed = seed))
  plugin_error <- ""
  plugin_asymptotic_se <- tryCatch({
    plugin_information <- information_matrix(
      fit$estimate, grid_n = config$plugin_grid_n
    )
    value <- sqrt(diag(safe_solve(plugin_information)))
    if (length(value) != length(theta_true) ||
        any(!is.finite(value)) || any(value <= 0))
      stop("invalid plug-in standard errors")
    names(value) <- names(theta_true)
    value
  }, error = function(e) {
    plugin_error <<- conditionMessage(e)
    rep(NA_real_, length(theta_true))
  })
  list(
    ok = TRUE, reason = "", seed = seed,
    rootn_error = sqrt(n) * (pack_theta(fit$estimate) - theta_true),
    plugin_asymptotic_se = plugin_asymptotic_se,
    plugin_error = plugin_error,
    H2 = fit$H2, grid_n = fit$grid_n,
    grid_step_over_h = fit$grid_step_over_h
  )
}

summarize_n <- function(replicates, n, truth_asymptotic_se, B) {
  valid <- vapply(replicates, function(result) isTRUE(result$ok), logical(1))
  if (sum(valid) < 2L)
    stop("Fewer than two valid replications at n=", n, ".")
  rootn <- do.call(rbind, lapply(
    replicates[valid], `[[`, "rootn_error"
  ))
  plugin_se <- do.call(rbind, lapply(
    replicates[valid], `[[`, "plugin_asymptotic_se"
  ))
  plugin_standardized <- rootn / plugin_se
  parameters <- colnames(rootn)
  rows <- lapply(seq_along(parameters), function(j) {
    plugin_ok <- is.finite(plugin_standardized[, j])
    data.frame(
      n = n, parameter = parameters[j],
      requested = B, valid = sum(valid), success_rate = mean(valid),
      z_bar = if (any(plugin_ok))
        mean(plugin_standardized[plugin_ok, j]) else NA_real_,
      empirical_sd = stats::sd(rootn[, j]) / sqrt(n),
      truth_sd = truth_asymptotic_se[j] / sqrt(n),
      plugin_sd = if (any(plugin_ok))
        mean(plugin_se[plugin_ok, j]) / sqrt(n) else NA_real_,
      plugin_coverage_95 = if (any(plugin_ok))
        mean(abs(plugin_standardized[plugin_ok, j]) <= stats::qnorm(0.975))
      else NA_real_,
      plugin_valid = sum(plugin_ok),
      stringsAsFactors = FALSE
    )
  })
  failures <- vapply(
    replicates[!valid], function(result) result$reason, character(1)
  )
  diagnostics <- data.frame(
    n = n, requested = B, valid = sum(valid),
    success_rate = mean(valid),
    plugin_valid = sum(vapply(
      replicates[valid], function(result)
        all(is.finite(result$plugin_asymptotic_se)), logical(1)
    )),
    plugin_success_rate = mean(vapply(
      replicates[valid], function(result)
        all(is.finite(result$plugin_asymptotic_se)), logical(1)
    )),
    failure_count = sum(!valid),
    first_failure = if (length(failures)) failures[1L] else "",
    first_plugin_failure = {
      plugin_failures <- vapply(
        replicates[valid], function(result) result$plugin_error, character(1)
      )
      plugin_failures <- plugin_failures[nzchar(plugin_failures)]
      if (length(plugin_failures)) plugin_failures[1L] else ""
    },
    stringsAsFactors = FALSE
  )
  list(summary = do.call(rbind, rows), diagnostics = diagnostics)
}

latex_parameter <- function(parameter) {
  stem <- sub("[0-9]+$", "", parameter)
  index <- sub("^[A-Za-z]+", "", parameter)
  symbol <- c(
    pi = "\\pi", mu = "\\mu", sigma = "\\sigma",
    alpha = "\\alpha", tau = "\\tau"
  )[[stem]]
  sprintf("$%s_%s$", symbol, index)
}

write_table3 <- function(summary, path, B, bandwidth_constant, beta) {
  format_signed <- function(x)
    sprintf("%+.2f", x)
  ns <- sort(unique(summary$n))
  parameters <- parameter_names(length(scenario3_truth$pi))
  lines <- c(
    sprintf(
      "%% Scenario 3 Table 3; B=%d; h_n=%g*n^{-%g}; generated by scenario3.R",
      B, bandwidth_constant, beta
    ),
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\begin{tabular}{l", paste(rep("ccccc", length(ns)), collapse = ""), "}"),
    "\\toprule",
    paste0(
      " & ",
      paste(sprintf("\\multicolumn{5}{c}{$n=%d$}", ns), collapse = " & "),
      " \\\\"
    ),
    paste(vapply(seq_along(ns), function(i)
      sprintf("\\cmidrule(lr){%d-%d}", 2L + 5L * (i - 1L), 6L + 5L * (i - 1L)),
      character(1)), collapse = " "),
    paste0(
      "$\\Theta$ & ",
      paste(rep(
        paste0(
          "$\\bar z$ & $\\mathrm{SD}_{\\rm emp}$",
          " & $\\mathrm{SD}_0$",
          " & $\\widehat{\\mathrm{SD}}_{\\rm plug}$",
          " & $\\mathrm{Cov}_{\\rm plug}$"
        ),
        length(ns)
      ), collapse = " & "),
      " \\\\"
    ),
    "\\midrule"
  )
  for (parameter in parameters) {
    cells <- character()
    for (n in ns) {
      row <- summary[
        summary$n == n & summary$parameter == parameter, , drop = FALSE
      ]
      if (nrow(row) != 1L)
        stop("Missing or duplicate table cell for ", parameter, ", n=", n)
      cells <- c(cells, sprintf(
        "%s & %.4f & %.4f & %.4f & %.3f",
        format_signed(row$z_bar), row$empirical_sd,
        row$truth_sd, row$plugin_sd, row$plugin_coverage_95
      ))
    }
    lines <- c(lines, paste0(
      latex_parameter(parameter), " & ", paste(cells, collapse = " & "), " \\\\"
    ))
  }
  lines <- c(
    lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"
  )
  writeLines(lines, path, useBytes = TRUE)
}

run_scenario3 <- function(
    B = 200L, n_values = c(500L, 1000L, 2000L), cores = 4L,
    seed = 20240717L, bandwidth_constant = 2, beta = 1 / 3,
    nstart = 10L, tolerance = 1e-5, max_iter = 100L,
    maxeval = 250L, kde_grid_n_max = 262144L,
    kde_points_per_bandwidth = 10, truth_grid_n = 200001L,
    plugin_grid_n = 20001L,
    active_weight_min = 0.01,
    out_dir = file.path(.scenario3_script_dir, "scenario3_results"),
    checkpoint_every = NULL, force = FALSE, quick = FALSE) {
  if (!requireNamespace("nloptr", quietly = TRUE))
    stop("R package 'nloptr' is required.")
  n_values <- as.integer(n_values)
  if (B < 1L || cores < 1L || any(!n_values %in% c(500L, 1000L, 2000L)))
    stop("Table 3 uses n=500,1000,2000; B and cores must be positive.")
  if (!(beta > 1 / (4 * min(pmin(scenario3_truth$alpha, 1))) &&
        beta < 0.5))
    stop("beta is outside the inference bandwidth window.")
  if (quick) {
    B <- 2L
    n_values <- 500L
    cores <- 1L
    nstart <- 3L
    max_iter <- 30L
    maxeval <- 100L
    truth_grid_n <- 20001L
    plugin_grid_n <- 5001L
  }
  if (is.null(checkpoint_every))
    checkpoint_every <- max(1L, 2L * min(B, cores))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(out_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  truth <- scenario3_truth
  theta_true <- pack_theta(truth)
  truth_information <- information_matrix(truth, grid_n = truth_grid_n)
  truth_covariance <- safe_solve(truth_information)
  truth_asymptotic_se <- sqrt(diag(truth_covariance))
  names(truth_asymptotic_se) <- names(theta_true)
  config <- list(
    version = "scenario3_table3_plugin_v2",
    B = B, n_values = n_values, seed = seed,
    bandwidth_constant = bandwidth_constant, beta = beta,
    nstart = nstart, tolerance = tolerance,
    max_iter = max_iter, maxeval = maxeval,
    kde_grid_n_max = kde_grid_n_max,
    kde_points_per_bandwidth = kde_points_per_bandwidth,
    truth_grid_n = truth_grid_n,
    plugin_grid_n = plugin_grid_n,
    active_weight_min = active_weight_min,
    truth = truth,
    standardization = paste(
      "(theta_hat-theta_true)/sqrt(diag(J(theta_hat)^-1)/n),",
      "equivalently sqrt(n)*(theta_hat-theta_true)/",
      "sqrt(diag(J(theta_hat)^-1))"
    ),
    variance = paste(
      "SD_0 reports J(theta_true)^-1/n for comparison only;",
      "standardization and coverage use J(theta_hat)^-1/n"
    )
  )
  summaries <- list()
  diagnostic_rows <- list()
  all_replicates <- list()

  for (n in n_values) {
    seed_base <- seed + n * 1000L
    metadata <- list(config = config, n = n, seed_base = seed_base)
    cache_file <- file.path(
      cache_dir, sprintf("scenario3_table3_plugin_v2_n%d_B%d.rds", n, B)
    )
    replicates <- vector("list", B)
    legacy_file <- file.path(
      cache_dir, sprintf("scenario3_table3_direct_v1_n%d_B%d.rds", n, B)
    )
    if (!force && !file.exists(cache_file) && file.exists(legacy_file)) {
      legacy <- readRDS(legacy_file)
      comparable <- c(
        "B", "n_values", "seed", "bandwidth_constant", "beta",
        "nstart", "tolerance", "max_iter", "maxeval",
        "kde_grid_n_max", "kde_points_per_bandwidth",
        "truth_grid_n", "plugin_grid_n", "active_weight_min", "truth"
      )
      compatible <- is.list(legacy$metadata$config) &&
        isTRUE(all.equal(
          legacy$metadata$config[comparable], config[comparable],
          tolerance = 0, check.attributes = FALSE
        )) &&
        identical(as.integer(legacy$metadata$n), as.integer(n)) &&
        identical(
          as.integer(legacy$metadata$seed_base), as.integer(seed_base)
        ) &&
        length(legacy$replicates) == B &&
        all(vapply(legacy$replicates, function(result)
          !is.null(result) && (
            !isTRUE(result$ok) ||
              length(result$plugin_asymptotic_se) == length(theta_true)
          ), logical(1)))
      if (compatible) {
        saveRDS(
          list(metadata = metadata, replicates = legacy$replicates),
          cache_file, compress = "gzip"
        )
        message("  migrated compatible completed v1 checkpoint: ",
                basename(legacy_file))
      }
    }
    if (!force && file.exists(cache_file)) {
      checkpoint <- readRDS(cache_file)
      if (!isTRUE(all.equal(
        checkpoint$metadata, metadata, tolerance = 0,
        check.attributes = FALSE
      ))) {
        stop(
          "Checkpoint metadata differ: ", cache_file,
          ". Use --force or a different --out directory."
        )
      }
      replicates <- checkpoint$replicates
      length(replicates) <- B
    }
    missing <- which(vapply(replicates, is.null, logical(1)))
    message(sprintf(
      "Scenario 3, n=%d: resumed %d/%d",
      n, B - length(missing), B
    ))

    if (length(missing)) {
      batches <- split(
        missing,
        ceiling(seq_along(missing) / as.integer(checkpoint_every))
      )
      cluster <- NULL
      if (cores > 1L && length(missing) > 1L) {
        cluster <- parallel::makeCluster(min(cores, length(missing)))
        scenario_file <- .scenario3_script_file
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
              n = n, truth = truth, theta_true = theta_true,
              config = config, seed_base = seed_base
            )
          } else {
            lapply(
              batch, one_replication,
              n = n, truth = truth, theta_true = theta_true,
              config = config, seed_base = seed_base
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

    result <- summarize_n(replicates, n, truth_asymptotic_se, B)
    summaries[[as.character(n)]] <- result$summary
    diagnostic_rows[[as.character(n)]] <- result$diagnostics
    all_replicates[[as.character(n)]] <- replicates
    message(sprintf(
      "  valid %d/%d; plug-in variance %d/%d; mean coverage %.3f",
      result$diagnostics$valid, B,
      result$diagnostics$plugin_valid, result$diagnostics$valid,
      mean(result$summary$plugin_coverage_95, na.rm = TRUE)
    ))
  }

  summary <- do.call(rbind, summaries)
  diagnostics <- do.call(rbind, diagnostic_rows)
  rownames(summary) <- rownames(diagnostics) <- NULL
  summary_file <- file.path(out_dir, "scenario3_summary.csv")
  diagnostic_file <- file.path(out_dir, "scenario3_diagnostics.csv")
  table_file <- file.path(out_dir, "scenario3_table3.tex")
  utils::write.csv(summary, summary_file, row.names = FALSE)
  utils::write.csv(diagnostics, diagnostic_file, row.names = FALSE)
  write_table3(summary, table_file, B, bandwidth_constant, beta)
  saveRDS(
    list(
      config = config,
      truth_information = truth_information,
      truth_covariance = truth_covariance,
      summary = summary, diagnostics = diagnostics,
      replicates = all_replicates
    ),
    file.path(out_dir, "scenario3_run_manifest.rds"),
    compress = "gzip"
  )
  message("Written:")
  message("  ", normalizePath(summary_file, winslash = "/"))
  message("  ", normalizePath(table_file, winslash = "/"))
  invisible(list(summary = summary, diagnostics = diagnostics))
}

parse_cli <- function(arguments) {
  config <- list(
    B = 200L, n_values = c(500L, 1000L, 2000L), cores = 4L,
    seed = 20240717L, bandwidth_constant = 2, beta = 1 / 3,
    nstart = 10L, tolerance = 1e-5, max_iter = 100L,
    maxeval = 250L, kde_grid_n_max = 262144L,
    kde_points_per_bandwidth = 10, truth_grid_n = 200001L,
    plugin_grid_n = 20001L,
    active_weight_min = 0.01,
    out_dir = file.path(.scenario3_script_dir, "scenario3_results"),
    checkpoint_every = NULL, force = FALSE, quick = FALSE
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
    else if (grepl("^--cores=", argument))
      config$cores <- as.integer(value_after(argument, "--cores="))
    else if (grepl("^--nstart=", argument))
      config$nstart <- as.integer(value_after(argument, "--nstart="))
    else if (grepl("^--seed=", argument))
      config$seed <- as.integer(value_after(argument, "--seed="))
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
  do.call(run_scenario3, cli)
}
