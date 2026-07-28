# Core routines for Scenario 1 in the PMHD-MCP paper.
#
# This file contains only the machinery used by scenario1.R:
#   * AEPD density, simulation, and numerical integration;
#   * the current block algorithm with the convex simplex weight update;
#   * exact-support and solver KKT diagnostics;
#   * regularization-path generation;
#   * Gaussian- and AEPD-likelihood baselines;
#   * density and parameter-error helpers.
#
# There is deliberately no pruning step. A zero weight can only be created by
# Euclidean projection in the convex weight subproblem. support_tol is a
# numerical reporting tolerance; it never modifies a fitted weight.

`%||%` <- function(x, y) if (is.null(x)) y else x

safe_scale <- function(x, fallback = 1) {
  candidates <- c(
    suppressWarnings(stats::IQR(x, na.rm = TRUE) / 1.349),
    suppressWarnings(stats::mad(x, constant = 1.4826, na.rm = TRUE)),
    suppressWarnings(stats::sd(x, na.rm = TRUE)),
    fallback
  )
  candidates <- candidates[is.finite(candidates) & candidates > 0]
  if (length(candidates)) candidates[1L] else 1
}

trapz_weights <- function(x) {
  if (length(x) < 2L || any(!is.finite(x)) || any(diff(x) <= 0))
    stop("The quadrature grid must contain at least two increasing values.")
  dx <- diff(x)
  w <- numeric(length(x))
  w[1L] <- dx[1L] / 2
  w[length(x)] <- dx[length(dx)] / 2
  if (length(x) > 2L)
    w[2L:(length(x) - 1L)] <- (dx[-length(dx)] + dx[-1L]) / 2
  w
}

normalize_grid_density <- function(y, w, floor = 0) {
  y <- pmax(as.numeric(y), floor)
  mass <- sum(y * w)
  if (!is.finite(mass) || mass <= 0)
    stop("A grid density has nonpositive or nonfinite mass.")
  y / mass
}

hellinger2_grid <- function(p, q, w, normalize = FALSE) {
  if (normalize) {
    p <- normalize_grid_density(p, w)
    q <- normalize_grid_density(q, w)
  }
  2 - 2 * sum(w * sqrt(pmax(p, 0) * pmax(q, 0)))
}

ise_grid <- function(p, q, w, normalize = FALSE) {
  if (normalize) {
    p <- normalize_grid_density(p, w)
    q <- normalize_grid_density(q, w)
  }
  sum(w * (p - q)^2)
}

# ---------------------------------------------------------------------------
# AEPD distribution
# ---------------------------------------------------------------------------

log_daepd <- function(x, mu, sigma, alpha, tau) {
  sigma <- pmax(sigma, 1e-12)
  alpha <- pmax(alpha, 1e-12)
  tau <- pmin(pmax(tau, 1e-10), 1 - 1e-10)
  side <- ifelse(x < mu, 1 - tau, tau)
  z <- (side * abs(x - mu) / sigma)^alpha
  log(alpha) + log(tau) + log1p(-tau) -
    lgamma(1 / alpha) - log(sigma) - z
}

daepd <- function(x, mu, sigma, alpha, tau) {
  sigma <- pmax(sigma, 1e-12)
  alpha <- pmax(alpha, 1e-12)
  tau <- pmin(pmax(tau, 1e-10), 1 - 1e-10)

  side <- tau + (1 - 2 * tau) * (x < mu)
  z <- (side * abs(x - mu) / sigma)^alpha
  cst <- alpha * tau * (1 - tau) /
    (gamma(1 / alpha) * sigma)

  pmax(cst * exp(-pmin(z, 745)), 1e-300)
}

raepd <- function(n, mu, sigma, alpha, tau) {
  right <- stats::runif(n) < (1 - tau)
  radius <- stats::rgamma(n, shape = 1 / alpha)^(1 / alpha) * sigma
  ifelse(right, mu + radius / tau, mu - radius / (1 - tau))
}

rmix_aepd <- function(n, pi, mu, sigma, alpha, tau) {
  K <- length(pi)
  if (!all(vapply(list(mu, sigma, alpha, tau), length, integer(1)) == K))
    stop("All mixture-parameter vectors must have the same length.")
  label <- sample.int(K, n, replace = TRUE, prob = pi)
  x <- numeric(n)
  for (k in seq_len(K)) {
    idx <- which(label == k)
    if (length(idx))
      x[idx] <- raepd(length(idx), mu[k], sigma[k], alpha[k], tau[k])
  }
  sample(x, length(x), replace = FALSE)
}

mixture_aepd_grid <- function(grid, pi, mu, sigma, alpha, tau) {
  G <- vapply(seq_along(pi), function(k)
    daepd(grid, mu[k], sigma[k], alpha[k], tau[k]),
    numeric(length(grid)))
  as.vector(G %*% pi)
}

aepd_component_moments <- function(mu, sigma, alpha, tau) {
  r1 <- gamma(2 / alpha) / gamma(1 / alpha)
  r2 <- gamma(3 / alpha) / gamma(1 / alpha)
  shift <- sigma * r1 * ((1 - tau) / tau - tau / (1 - tau))
  raw2 <- sigma^2 * r2 *
    ((1 - tau) / tau^2 + tau / (1 - tau)^2)
  list(mean = mu + shift, sd = sqrt(pmax(raw2 - shift^2, 0)))
}

# R's bandwidth convention scales the Epanechnikov kernel to unit variance.
epanechnikov_kernel <- function(u) {
  a <- sqrt(5)
  out <- numeric(length(u))
  keep <- abs(u) <= a
  out[keep] <- 3 / (4 * a) * (1 - u[keep]^2 / 5)
  out
}

make_kde_grid <- function(x, bw_adjust = 0.80, grid_n = 512L) {
  x <- as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x)))
    stop("x must contain at least two finite observations.")
  h <- bw_adjust * stats::bw.nrd0(x)
  if (!is.finite(h) || h <= 0)
    h <- bw_adjust * 1.06 * safe_scale(x) * length(x)^(-1 / 5)
  support <- sqrt(5) * h
  grid <- seq(min(x) - support, max(x) + support, length.out = grid_n)
  u <- outer(grid, x, "-") / h
  phat <- rowMeans(matrix(epanechnikov_kernel(u), nrow = length(grid))) / h
  w <- trapz_weights(grid)
  phat <- normalize_grid_density(phat, w)
  list(
    grid = grid, w = w, phat = phat, sqrt_phat = sqrt(phat),
    bandwidth = h, bw_adjust = bw_adjust, kernel = "epanechnikov",
    grid_n = as.integer(grid_n), grid_range = range(grid)
  )
}

# ---------------------------------------------------------------------------
# MCP and the convex weight subproblem
# ---------------------------------------------------------------------------

mcp_value <- function(pi, lambda, gamma) {
  ifelse(
    pi <= gamma * lambda,
    lambda * pi - pi^2 / (2 * gamma),
    gamma * lambda^2 / 2
  )
}

mcp_derivative <- function(pi, lambda, gamma) {
  pmax(lambda - pi / gamma, 0)
}

project_simplex_mass <- function(v, mass = 1) {
  if (!length(v) || any(!is.finite(v)) || !is.finite(mass) || mass <= 0)
    stop("Invalid simplex-projection input.")
  u <- sort(v, decreasing = TRUE)
  cs <- cumsum(u)
  idx <- which(u - (cs - mass) / seq_along(u) > 0)
  if (!length(idx))
    stop("The simplex projection failed.")
  rho <- max(idx)
  theta <- (cs[rho] - mass) / rho
  p <- pmax(v - theta, 0)
  positive <- which(p > 0)
  drift <- mass - sum(p)
  if (length(positive))
    p[positive[which.max(p[positive])]] <-
      p[positive[which.max(p[positive])]] + drift
  p[p < 0 & p > -100 * .Machine$double.eps] <- 0
  p
}

project_simplex <- function(v, lower = 0) {
  lower <- rep_len(lower, length(v))
  if (any(!is.finite(lower)) || any(lower < 0) || sum(lower) >= 1)
    stop("lower must be nonnegative and sum to less than one.")
  lower + project_simplex_mass(v - lower, 1 - sum(lower))
}

component_matrix <- function(grid, pars) {
  vapply(seq_along(pars$pi), function(k)
    daepd(grid, pars$mu[k], pars$sigma[k], pars$alpha[k], pars$tau[k]),
    numeric(length(grid)))
}

hellinger2_weights <- function(pi, sqrt_phat, G, w) {
  f <- as.vector(G %*% pi)
  2 - 2 * sum(w * sqrt_phat * sqrt(pmax(f, 0)))
}

G_vector <- function(pi, sqrt_phat, G, w, floor = 1e-300) {
  f <- as.vector(G %*% pi)
  as.vector(crossprod(G, w * sqrt_phat / sqrt(pmax(f, floor))))
}

solve_weight_convex <- function(pi0, sqrt_phat, G, w, xi,
                                lower = 0, max_iter = 300L,
                                pg_tol = 1e-6, verbose = FALSE) {
  K <- ncol(G)
  if (length(pi0) != K || length(xi) != K)
    stop("pi0 and xi must match the number of component columns.")
  p <- project_simplex(pi0, lower)
  objective <- function(q)
    hellinger2_weights(q, sqrt_phat, G, w) + sum(xi * q)
  gradient <- function(q) -G_vector(q, sqrt_phat, G, w) + xi
  fp <- objective(p)
  g <- gradient(p)
  step_size <- 1 / max(1, max(abs(g)))
  pg <- max(abs(p - project_simplex(p - g, lower)))
  reason <- "maximum iterations reached"
  converged <- pg <= pg_tol
  iterations <- 0L

  for (iteration in seq_len(max_iter)) {
    if (converged) break
    iterations <- iteration
    accepted <- FALSE
    for (line_search in seq_len(60L)) {
      candidate <- project_simplex(p - step_size * g, lower)
      direction <- candidate - p
      fc <- objective(candidate)
      upper <- fp + sum(g * direction) +
        sum(direction^2) / (2 * step_size)
      if (is.finite(fc) && fc <= upper + 1e-14) {
        accepted <- TRUE
        break
      }
      step_size <- step_size / 2
    }
    if (!accepted) {
      reason <- "backtracking failed"
      break
    }
    p <- candidate
    fp <- fc
    g <- gradient(p)
    pg <- max(abs(p - project_simplex(p - g, lower)))
    step_size <- min(step_size * 2, 1e6)
    if (verbose)
      message(sprintf("weight iteration %d: objective %.10g, PG %.3e",
                      iteration, fp, pg))
    if (pg <= pg_tol) {
      converged <- TRUE
      reason <- "projected-gradient tolerance satisfied"
    }
  }

  list(
    pi = p, objective = fp, converged = converged,
    iterations = iterations, reason = reason,
    projected_gradient_residual = pg,
    exact_zero_count = sum(p == 0)
  )
}

solve_weight_lla <- function(pi0, sqrt_phat, G, w, lambda, gamma,
                             max_lla_iter = 10L,
                             convex_max_iter = 300L,
                             pg_tol = 1e-6, kkt_tol = 1e-5,
                             support_tol = 1e-8) {
  # Repeating the LLA weight update with theta fixed is a sequence of valid
  # MM steps: the unchanged theta is an admissible (equal-affinity) theta
  # update. This cheaply brings the exact MCP weight KKT conditions to the
  # requested tolerance before another expensive Nelder--Mead sweep.
  p <- project_simplex(pi0)
  pg_trace <- numeric(0)
  convex_iterations <- integer(0)
  final_kkt <- NULL
  final_step <- NULL
  reason <- "maximum LLA weight cycles reached"

  for (lla_iteration in seq_len(max_lla_iter)) {
    final_step <- solve_weight_convex(
      p, sqrt_phat, G, w, mcp_derivative(p, lambda, gamma),
      max_iter = convex_max_iter, pg_tol = pg_tol
    )
    p <- final_step$pi
    final_kkt <- weight_kkt_check(
      p, sqrt_phat, G, w, lambda, gamma,
      tolerance = kkt_tol, support_tol = support_tol
    )
    pg_trace <- c(pg_trace, final_kkt$projected_gradient_residual)
    convex_iterations <- c(convex_iterations, final_step$iterations)
    if (isTRUE(final_kkt$certified)) {
      reason <- "exact MCP weight KKT tolerance satisfied"
      break
    }
  }

  list(
    pi = p,
    converged = isTRUE(final_kkt$certified),
    reason = reason,
    lla_iterations = length(pg_trace),
    convex_iterations = convex_iterations,
    projected_gradient_residual =
      final_kkt$projected_gradient_residual,
    exact_zero_count = sum(p == 0),
    kkt = final_kkt,
    pg_trace = pg_trace,
    last_convex_step = final_step
  )
}

weight_kkt_check <- function(pi, sqrt_phat, G, w, lambda, gamma,
                             tolerance = 1e-5,
                             support_tol = 1e-8) {
  Gv <- G_vector(pi, sqrt_phat, G, w)
  rho <- mcp_derivative(pi, lambda, gamma)
  active <- which(pi > 0)
  zero <- which(pi == 0)
  near <- which(pi > 0 & pi <= support_tol)
  if (!length(active))
    stop("A simplex point must have at least one active coordinate.")

  active_gradient <- rho[active] - Gv[active]
  nu <- mean(active_gradient)
  active_residual <- max(abs(active_gradient - nu))
  zero_slack <- if (length(zero))
    lambda - Gv[zero] - nu else numeric(0)
  zero_residual <- if (length(zero))
    max(pmax(-zero_slack, 0)) else 0
  full_gradient <- -Gv + rho
  pg_residual <- max(abs(pi - project_simplex(pi - full_gradient)))
  simplex_residual <- abs(sum(pi) - 1)
  nonnegative_residual <- max(pmax(-pi, 0))
  certified <- all(is.finite(c(
    active_residual, zero_residual, pg_residual,
    simplex_residual, nonnegative_residual
  ))) &&
    active_residual <= tolerance &&
    zero_residual <= tolerance &&
    pg_residual <= tolerance &&
    simplex_residual <= tolerance &&
    nonnegative_residual <= tolerance

  list(
    certified = certified, nu = nu, active = active, zero = zero,
    near_zero_positive = near, active_residual = active_residual,
    zero_residual = zero_residual, zero_slack = zero_slack,
    projected_gradient_residual = pg_residual,
    simplex_residual = simplex_residual,
    nonnegative_residual = nonnegative_residual,
    G = Gv
  )
}

# ---------------------------------------------------------------------------
# Initial values and single-component Hellinger updates
# ---------------------------------------------------------------------------

make_parameter_bounds <- function(x, alpha_lower = 1, alpha_upper = 3,
                                  tau_lower = 0.01, tau_upper = 0.99,
                                  grid_spacing = NULL) {
  sx <- safe_scale(x)
  rx <- diff(range(x))
  if (!is.finite(rx) || rx <= 0) rx <- 4 * sx
  sigma_lower <- max(1e-3 * sx, 1e-5)
  if (!is.null(grid_spacing)) {
    if (!is.finite(grid_spacing) || grid_spacing <= 0)
      stop("grid_spacing must be positive when supplied.")
    # A narrower component is not resolved by the common quadrature grid and
    # can acquire spurious mass when its mode falls on a grid point.
    sigma_lower <- max(sigma_lower, 2 * grid_spacing)
  }
  list(
    mu = c(min(x) - sx, max(x) + sx),
    log_sigma = log(c(sigma_lower, max(3 * rx, sx))),
    log_alpha = log(c(alpha_lower, alpha_upper)),
    logit_tau = stats::qlogis(c(tau_lower, tau_upper)),
    alpha_lower = alpha_lower, alpha_upper = alpha_upper,
    tau_lower = tau_lower, tau_upper = tau_upper
  )
}

clip <- function(x, bounds, margin = 1e-10) {
  pmin(pmax(x, bounds[1L] + margin * diff(bounds)),
       bounds[2L] - margin * diff(bounds))
}

component_to_free <- function(par, bounds) {
  transformed <- c(
    par["mu"], log(par["sigma"]), log(par["alpha"]),
    stats::qlogis(par["tau"])
  )
  box <- rbind(
    bounds$mu, bounds$log_sigma, bounds$log_alpha, bounds$logit_tau
  )
  unit <- vapply(seq_len(4L), function(j)
    (clip(transformed[j], box[j, ]) - box[j, 1L]) /
      diff(box[j, ]), numeric(1))
  stats::qlogis(pmin(pmax(unit, 1e-8), 1 - 1e-8))
}

free_to_component <- function(z, bounds) {
  box <- rbind(
    bounds$mu, bounds$log_sigma, bounds$log_alpha, bounds$logit_tau
  )
  transformed <- box[, 1L] +
    (box[, 2L] - box[, 1L]) * stats::plogis(z)
  c(
    mu = transformed[1L],
    sigma = exp(transformed[2L]),
    alpha = exp(transformed[3L]),
    tau = stats::plogis(transformed[4L])
  )
}

sort_parameters <- function(pars) {
  ord <- order(pars$mu, pars$alpha, pars$tau)
  lapply(pars, function(v) v[ord])
}

sanitize_initial <- function(init, K, bounds) {
  required <- c("pi", "mu", "sigma", "alpha", "tau")
  if (!all(required %in% names(init)) ||
      any(vapply(init[required], length, integer(1)) != K))
    stop("The initial value is incomplete or has the wrong order.")
  init$pi <- project_simplex(pmax(init$pi, 0))
  init$mu <- clip(init$mu, bounds$mu)
  init$sigma <- exp(clip(log(pmax(init$sigma, 1e-12)),
                           bounds$log_sigma))
  init$alpha <- exp(clip(log(pmax(init$alpha, 1e-12)),
                           bounds$log_alpha))
  init$tau <- stats::plogis(clip(stats::qlogis(
    pmin(pmax(init$tau, 1e-10), 1 - 1e-10)), bounds$logit_tau))
  sort_parameters(init[required])
}

initial_aepd_mixture <- function(x, K, bounds, random = FALSE) {
  sx <- safe_scale(x)
  km <- tryCatch(
    stats::kmeans(matrix(x, ncol = 1L), centers = K,
                  nstart = if (random) 1L else 20L),
    error = function(e) NULL
  )
  if (is.null(km)) {
    centers <- as.numeric(stats::quantile(
      x, (seq_len(K) - 0.5) / K, names = FALSE, type = 8
    ))
    cluster <- max.col(-abs(outer(x, centers, "-")), ties.method = "first")
  } else {
    ord <- order(km$centers[, 1L])
    centers <- km$centers[ord, 1L]
    cluster <- match(km$cluster, ord)
  }
  pi <- pmax(tabulate(cluster, nbins = K) / length(x), 1e-6)
  mu <- vapply(seq_len(K), function(k) {
    z <- x[cluster == k]
    if (length(z)) stats::median(z) else centers[k]
  }, numeric(1))
  sigma <- vapply(seq_len(K), function(k) {
    z <- x[cluster == k]
    max(safe_scale(z, sx / max(K, 2)) / sqrt(2), 0.03 * sx)
  }, numeric(1))
  tau <- vapply(seq_len(K), function(k) {
    z <- x[cluster == k]
    if (length(z) < 5L) return(0.5)
    left <- safe_scale(mu[k] - z[z <= mu[k]], sigma[k])
    right <- safe_scale(z[z > mu[k]] - mu[k], sigma[k])
    pmin(pmax(left / (left + right), 0.15), 0.85)
  }, numeric(1))
  init <- list(
    pi = pi / sum(pi), mu = mu, sigma = sigma,
    alpha = rep(2, K), tau = tau
  )
  if (random) {
    init$pi <- stats::rgamma(K, shape = 2 + 10 * init$pi)
    init$pi <- init$pi / sum(init$pi)
    init$mu <- init$mu + stats::rnorm(K, sd = 0.10 * sx)
    init$sigma <- init$sigma * exp(stats::rnorm(K, sd = 0.25))
    init$alpha <- stats::runif(
      K, bounds$alpha_lower, bounds$alpha_upper
    )
    init$tau <- stats::runif(
      K, max(bounds$tau_lower, 0.15), min(bounds$tau_upper, 0.85)
    )
  }
  sanitize_initial(init, K, bounds)
}

optimize_component_affinity <- function(start, target, component_function,
                                        w, bounds,
                                        stage1_maxit = 120L,
                                        stage2_maxit = 240L) {
  target <- pmax(target, 0)
  affinity <- function(par)
    2 * sum(w * sqrt(target * pmax(component_function(par), 0)))
  objective <- function(z) -affinity(free_to_component(z, bounds))
  z0 <- component_to_free(start, bounds)
  candidates <- list(start)
  first <- tryCatch(
    stats::optim(
      z0, objective, method = "Nelder-Mead",
      control = list(maxit = stage1_maxit, reltol = 1e-6)
    ),
    error = function(e) NULL
  )
  if (!is.null(first) && all(is.finite(first$par))) {
    p1 <- free_to_component(first$par, bounds)
    candidates[[length(candidates) + 1L]] <- p1
    second <- tryCatch(
      stats::optim(
        first$par, objective, method = "Nelder-Mead",
        control = list(maxit = stage2_maxit, reltol = 1e-9)
      ),
      error = function(e) NULL
    )
    if (!is.null(second) && all(is.finite(second$par)))
      candidates[[length(candidates) + 1L]] <-
        free_to_component(second$par, bounds)
  }
  values <- vapply(candidates, affinity, numeric(1))
  best <- which.max(values)
  list(par = candidates[[best]], affinity = values[best])
}

parameter_change <- function(old, new, scale_mu) {
  max(c(
    abs(old$pi - new$pi),
    abs(old$mu - new$mu) / scale_mu,
    abs(log(old$sigma) - log(new$sigma)),
    abs(log(old$alpha) - log(new$alpha)),
    abs(stats::qlogis(old$tau) - stats::qlogis(new$tau))
  ))
}

near_coincident_components <- function(pars, support_tol, scale_mu,
                                       tolerance = 1e-4) {
  active <- which(pars$pi > support_tol)
  if (length(active) < 2L)
    return(list(flag = FALSE, minimum_distance = Inf, pair = integer(0)))
  pairs <- utils::combn(active, 2L)
  distance <- apply(pairs, 2L, function(idx) sqrt(
    ((pars$mu[idx[1L]] - pars$mu[idx[2L]]) / scale_mu)^2 +
      ((pars$alpha[idx[1L]] - pars$alpha[idx[2L]]) /
         max(1, diff(range(pars$alpha))))^2
  ))
  j <- which.min(distance)
  list(
    flag = distance[j] <= tolerance,
    minimum_distance = distance[j], pair = pairs[, j]
  )
}

# ---------------------------------------------------------------------------
# PMHD-MCP block algorithm
# ---------------------------------------------------------------------------

pmhd_mcp_fit <- function(x, K, lambda, kde, init = NULL, gamma = 3,
                         alpha_lower = 1, alpha_upper = 3,
                         tau_lower = 0.01, tau_upper = 0.99,
                         max_outer = 40L,
                         theta_stage1 = 30L, theta_stage2 = 60L,
                         weight_max_iter = 300L,
                         weight_lla_max_iter = 10L,
                         parameter_tol = 1e-2, objective_tol = 2e-6,
                         weight_pg_tol = 1e-6, kkt_tol = 1e-5,
                         support_tol = 1e-8,
                         descent_tol = 1e-8,
                         coincidence_tol = 1e-4,
                         verbose = FALSE) {
  if (!is.finite(lambda) || lambda < 0 || !is.finite(gamma) || gamma <= 1)
    stop("lambda must be nonnegative and gamma must exceed one.")
  bounds <- make_parameter_bounds(
    x, alpha_lower, alpha_upper, tau_lower, tau_upper,
    grid_spacing = mean(diff(kde$grid))
  )
  if (is.null(init))
    init <- initial_aepd_mixture(x, K, bounds, random = FALSE)
  pars <- sanitize_initial(init, K, bounds)
  grid <- kde$grid
  w <- kde$w
  sqrt_phat <- kde$sqrt_phat
  phat <- kde$phat
  mu_scale <- safe_scale(x)

  G <- component_matrix(grid, pars)
  H2 <- hellinger2_weights(pars$pi, sqrt_phat, G, w)
  Q <- H2 + sum(mcp_value(pars$pi, lambda, gamma))
  objective_trace <- Q
  inner_pg_trace <- numeric(0)
  inner_lla_trace <- integer(0)
  accepted_theta_trace <- integer(0)
  last_change <- Inf
  last_objective_change <- Inf
  reason <- "maximum outer iterations reached"
  descent_ok <- TRUE
  outer_iterations <- 0L

  for (iteration in seq_len(max_outer)) {
    outer_iterations <- iteration
    old <- pars
    old_G <- G
    old_Q <- Q
    f_old <- pmax(as.vector(old_G %*% old$pi), 1e-300)
    responsibility <- sweep(
      sweep(old_G, 2L, old$pi, "*"), 1L, f_old, "/"
    )

    candidate <- old
    accepted_count <- 0L
    for (k in which(old$pi > 0)) {
      target <- phat * responsibility[, k]
      if (sum(target * w) <= 1e-14) next
      start <- c(
        mu = old$mu[k], sigma = old$sigma[k],
        alpha = old$alpha[k], tau = old$tau[k]
      )
      component_function <- function(par)
        daepd(grid, par["mu"], par["sigma"], par["alpha"], par["tau"])
      update <- optimize_component_affinity(
        start, target, component_function, w, bounds,
        stage1_maxit = theta_stage1, stage2_maxit = theta_stage2
      )
      old_affinity <- 2 * sum(w * sqrt(target * old_G[, k]))
      if (is.finite(update$affinity) &&
          update$affinity + 1e-12 >= old_affinity) {
        candidate$mu[k] <- update$par["mu"]
        candidate$sigma[k] <- update$par["sigma"]
        candidate$alpha[k] <- update$par["alpha"]
        candidate$tau[k] <- update$par["tau"]
        accepted_count <- accepted_count + 1L
      }
    }
    G_candidate <- component_matrix(grid, candidate)
    weight_update <- solve_weight_lla(
      old$pi, sqrt_phat, G_candidate, w, lambda, gamma,
      max_lla_iter = weight_lla_max_iter,
      convex_max_iter = weight_max_iter,
      pg_tol = weight_pg_tol, kkt_tol = kkt_tol,
      support_tol = support_tol
    )
    candidate$pi <- weight_update$pi
    # Sorting only permutes component labels. Reorder the already evaluated
    # component matrix instead of evaluating every density a second time.
    ord <- order(candidate$mu, candidate$alpha, candidate$tau)
    candidate <- lapply(candidate, function(v) v[ord])
    G_candidate <- G_candidate[, ord, drop = FALSE]
    H2_candidate <- hellinger2_weights(
      candidate$pi, sqrt_phat, G_candidate, w
    )
    Q_candidate <- H2_candidate +
      sum(mcp_value(candidate$pi, lambda, gamma))

    if (!is.finite(Q_candidate) || Q_candidate > old_Q + descent_tol) {
      descent_ok <- FALSE
      reason <- "the discretized objective increased"
      break
    }

    pars <- candidate
    G <- G_candidate
    H2 <- H2_candidate
    Q <- Q_candidate
    last_change <- parameter_change(old, pars, mu_scale)
    last_objective_change <- abs(old_Q - Q)
    objective_trace <- c(objective_trace, Q)
    inner_pg_trace <- c(
      inner_pg_trace, weight_update$projected_gradient_residual
    )
    inner_lla_trace <- c(
      inner_lla_trace, weight_update$lla_iterations
    )
    accepted_theta_trace <- c(accepted_theta_trace, accepted_count)
    kkt_now <- weight_kkt_check(
      pars$pi, sqrt_phat, G, w, lambda, gamma,
      tolerance = kkt_tol, support_tol = support_tol
    )

    if (verbose)
      message(sprintf(
        "outer %d: Q %.9g, K %d, change %.2e, KKT %.2e",
        iteration, Q, sum(pars$pi > support_tol), last_change,
        max(kkt_now$active_residual, kkt_now$zero_residual,
            kkt_now$projected_gradient_residual)
      ))

    if (last_change <= parameter_tol &&
        last_objective_change <= objective_tol &&
        kkt_now$certified) {
      reason <- "all numerical tolerances satisfied"
      break
    }
  }

  kkt <- weight_kkt_check(
    pars$pi, sqrt_phat, G, w, lambda, gamma,
    tolerance = kkt_tol, support_tol = support_tol
  )
  coincidence <- near_coincident_components(
    pars, support_tol, mu_scale, coincidence_tol
  )
  converged <- identical(reason, "all numerical tolerances satisfied") &&
    descent_ok && kkt$certified
  admissible <- converged && !coincidence$flag

  c(pars, list(
    K_working = K,
    K_exact = sum(pars$pi > 0),
    Khat = sum(pars$pi > support_tol),
    exact_zero_count = sum(pars$pi == 0),
    near_zero_positive = which(pars$pi > 0 & pars$pi <= support_tol),
    lambda = lambda, gamma = gamma, H2 = H2, Q = Q,
    converged = converged, admissible = admissible, reason = reason,
    descent_ok = descent_ok, outer_iterations = outer_iterations,
    parameter_change = last_change,
    objective_change = last_objective_change,
    objective_trace = objective_trace,
    inner_projected_gradient_trace = inner_pg_trace,
    inner_lla_iterations_trace = inner_lla_trace,
    accepted_theta_trace = accepted_theta_trace,
    kkt = kkt, coincidence = coincidence,
    support_tol = support_tol, bounds = bounds,
    bandwidth = kde$bandwidth, grid_range = kde$grid_range,
    grid_n = kde$grid_n
  ))
}

active_part <- function(fit, support_tol = fit$support_tol %||% 1e-8) {
  keep <- which(fit$pi > support_tol)
  if (!length(keep))
    stop("The fitted model has no reported active component.")
  out <- lapply(fit[c("pi", "mu", "sigma", "alpha", "tau")],
                function(v) v[keep])
  out$pi <- out$pi / sum(out$pi)
  out$Khat <- length(keep)
  out
}

make_start_list <- function(x, K, nstart, bounds, warm = NULL, seed = 1L) {
  set.seed(seed)
  starts <- vector("list", nstart)
  position <- 1L
  if (!is.null(warm)) {
    warm_list <- if (!is.null(warm$pi)) list(warm) else warm
    # Preserve pathwise branches without giving up fresh exploration at the
    # new lambda. At most half of the starts are warm starts.
    warm_limit <- max(1L, floor(nstart / 2L))
    warm_used <- 0L
    for (item in warm_list) {
      if (position > nstart || warm_used >= warm_limit) break
      candidate <- tryCatch(
        sanitize_initial(item[c("pi", "mu", "sigma", "alpha", "tau")],
                         K, bounds),
        error = function(e) NULL
      )
      if (!is.null(candidate)) {
        starts[[position]] <- candidate
        position <- position + 1L
        warm_used <- warm_used + 1L
      }
    }
  }
  # At the first lambda include one deterministic full-support start. Along
  # the path, a warm start already fills that role, so reserve the remaining
  # slots for genuinely fresh exploration.
  if (position <= nstart && is.null(warm)) {
    starts[[position]] <- initial_aepd_mixture(x, K, bounds, random = FALSE)
    position <- position + 1L
  }
  fresh_index <- 0L
  while (position <= nstart) {
    fresh_index <- fresh_index + 1L
    candidate <- initial_aepd_mixture(x, K, bounds, random = TRUE)
    # Cycle through boundary supports 1,...,K. These are initial values, not
    # pruning: the convex simplex update can reactivate every zero coordinate.
    support_size <- 1L + ((fresh_index - 1L) %% K)
    if (support_size < K) {
      keep <- sample.int(K, support_size)
      candidate$pi[-keep] <- 0
      candidate$pi <- candidate$pi / sum(candidate$pi)
    }
    starts[[position]] <- candidate
    position <- position + 1L
  }
  starts
}

fit_pmhd_at_lambda <- function(x, K, lambda, kde, nstart, gamma,
                               warm = NULL, seed = 1L,
                               fit_control = list()) {
  bounds <- make_parameter_bounds(
    x,
    fit_control$alpha_lower %||% 1,
    fit_control$alpha_upper %||% 3,
    fit_control$tau_lower %||% 0.01,
    fit_control$tau_upper %||% 0.99,
    grid_spacing = mean(diff(kde$grid))
  )
  starts <- make_start_list(x, K, nstart, bounds, warm, seed)
  fits <- lapply(seq_along(starts), function(s) {
    args <- c(list(
      x = x, K = K, lambda = lambda, kde = kde,
      init = starts[[s]], gamma = gamma
    ), fit_control)
    tryCatch(do.call(pmhd_mcp_fit, args), error = function(e)
      structure(list(error = conditionMessage(e)), class = "pmhd_error"))
  })
  good <- which(vapply(
    fits, function(z) is.list(z) && isTRUE(z$admissible) &&
      is.finite(z$Q), logical(1)
  ))

  # A single numerical retry is allowed; it changes no statistical rule.
  if (!length(good)) {
    finite <- which(vapply(
      fits, function(z) is.list(z) && !inherits(z, "pmhd_error") &&
        is.finite(z$Q %||% NA_real_), logical(1)
    ))
    if (length(finite)) {
      best <- finite[which.min(vapply(
        fits[finite], function(z) z$Q, numeric(1)
      ))]
      retry_control <- fit_control
      retry_control$max_outer <-
        (fit_control$max_outer %||% 40L) + 20L
      retry <- tryCatch(do.call(pmhd_mcp_fit, c(list(
        x = x, K = K, lambda = lambda, kde = kde,
        init = fits[[best]], gamma = gamma
      ), retry_control)), error = function(e) NULL)
      if (!is.null(retry)) {
        fits[[length(fits) + 1L]] <- retry
        if (isTRUE(retry$admissible)) good <- length(fits)
      }
    }
  }

  if (!length(good)) {
    return(list(
      success = FALSE, fit = NULL, fits = fits, lambda = lambda,
      reason = "no random start passed convergence, KKT, and separation checks"
    ))
  }
  Q <- vapply(fits[good], function(z) z$Q, numeric(1))
  best <- good[which.min(Q)]
  list(
    success = TRUE, fit = fits[[best]], fits = fits,
    lambda = lambda, retained_start = best,
    n_admissible = length(good), reason = "success"
  )
}

monotone_order_path <- function(K_raw) {
  out <- rep(NA_integer_, length(K_raw))
  current <- Inf
  for (i in seq_along(K_raw)) {
    if (is.finite(K_raw[i])) {
      current <- min(current, K_raw[i])
      out[i] <- as.integer(current)
    }
  }
  out
}

compress_h_path <- function(path) {
  K_levels <- unique(path$K_monotone[is.finite(path$K_monotone)])
  positive_lambda <- path$lambda[is.finite(path$lambda) & path$lambda > 0]
  lambda_floor <- if (length(positive_lambda))
    min(positive_lambda) / 2 else 1e-12
  rows <- lapply(K_levels, function(K) {
    eligible <- which(
      path$K_monotone == K & path$K_raw == K & is.finite(path$H2)
    )
    if (!length(eligible)) return(NULL)
    best <- eligible[which.min(path$H2[eligible])]
    log_lambda <- log(pmax(path$lambda[eligible], lambda_floor))
    data.frame(
      K = K, H2_star = path$H2[best],
      lambda_star = path$lambda[best],
      path_row = best, key = path$key[best],
      plateau_n = length(eligible),
      log_lambda_range = max(log_lambda) - min(log_lambda),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows))
    return(data.frame(
      K = integer(0), H2_star = numeric(0), lambda_star = numeric(0),
      path_row = integer(0), key = character(0),
      plateau_n = integer(0), log_lambda_range = numeric(0)
    ))
  do.call(rbind, rows)
}

default_lambda_grid <- function() {
  exp(seq(log(0.0015), log(0.30), length.out = 15L))
}

fit_pmhd_path <- function(x, K, lambda_grid = default_lambda_grid(),
                          nstart = 20L, gamma = 3,
                          bw_adjust = 0.80, grid_n = 512L,
                          seed = 1L, support_tol = 1e-8,
                          max_grid_expansions = 3L,
                          adaptive_nstart = min(5L, nstart),
                          fit_control = list()) {
  lambda_grid <- sort(unique(as.numeric(lambda_grid)))
  if (!length(lambda_grid) || any(!is.finite(lambda_grid)) ||
      any(lambda_grid < 0))
    stop("lambda_grid must contain finite nonnegative values.")
  if (!is.finite(nstart) || nstart < 1L)
    stop("nstart must be a positive integer.")
  if (!is.finite(max_grid_expansions) || max_grid_expansions < 0L)
    stop("max_grid_expansions must be a nonnegative integer.")
  if (!is.finite(adaptive_nstart) || adaptive_nstart < 1L)
    stop("adaptive_nstart must be a positive integer.")
  nstart <- as.integer(nstart)
  max_grid_expansions <- as.integer(max_grid_expansions)
  kde <- make_kde_grid(x, bw_adjust, grid_n)
  fit_control$support_tol <- fit_control$support_tol %||% support_tol
  records <- list()

  key_of <- function(lambda) sprintf("%.17g", lambda)
  adaptive_nstart <- max(
    1L, min(as.integer(adaptive_nstart), as.integer(nstart))
  )
  fit_missing_lambda <- function(lambda, starts_at_lambda = nstart) {
    key <- key_of(lambda)
    if (!is.null(records[[key]])) return(invisible(NULL))
    existing <- Filter(function(z) isTRUE(z$success), records)
    warm <- NULL
    if (length(existing)) {
      distance <- vapply(existing, function(z)
        abs(log1p(z$lambda) - log1p(lambda)), numeric(1))
      nearest <- existing[[which.min(distance)]]
      warm <- Filter(function(z)
        is.list(z) && !inherits(z, "pmhd_error") &&
          is.finite(z$Q %||% NA_real_), nearest$fits)
      if (!length(warm)) warm <- nearest$fit
    }
    record <- fit_pmhd_at_lambda(
      x, K, lambda, kde, starts_at_lambda, gamma, warm,
      seed = seed + as.integer(1e6 * lambda) + length(records) * 1009L,
      fit_control = fit_control
    )
    records[[key]] <<- record
    invisible(NULL)
  }
  build_path <- function() {
    lambdas <- sort(vapply(records, function(z) z$lambda, numeric(1)))
    rows <- lapply(lambdas, function(lambda) {
      key <- key_of(lambda)
      z <- records[[key]]
      if (!isTRUE(z$success)) {
        return(data.frame(
          key = key, lambda = lambda, success = FALSE,
          K_exact = NA_integer_, K_raw = NA_integer_,
          exact_zeros = NA_integer_, near_zero_count = NA_integer_,
          H2 = NA_real_, Q = NA_real_, active_kkt = NA_real_,
          zero_kkt = NA_real_, pg_kkt = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      f <- z$fit
      data.frame(
        key = key, lambda = lambda, success = TRUE,
        K_exact = f$K_exact, K_raw = f$Khat,
        exact_zeros = f$exact_zero_count,
        near_zero_count = length(f$near_zero_positive),
        H2 = f$H2, Q = f$Q,
        active_kkt = f$kkt$active_residual,
        zero_kkt = f$kkt$zero_residual,
        pg_kkt = f$kkt$projected_gradient_residual,
        stringsAsFactors = FALSE
      )
    })
    path <- do.call(rbind, rows)
    path$K_monotone <- monotone_order_path(path$K_raw)
    path
  }

  for (lambda in lambda_grid) fit_missing_lambda(lambda)
  path <- build_path()
  compressed <- compress_h_path(path)

  expansion <- 0L
  while (expansion < max_grid_expansions) {
    successful <- path[path$success, , drop = FALSE]
    new_lambda <- numeric(0)
    if (nrow(successful) >= 2L) {
      transitions <- which(diff(successful$K_monotone) != 0L)
      for (j in transitions) {
        lo <- successful$lambda[j]
        hi <- successful$lambda[j + 1L]
        order_gap <- abs(
          successful$K_monotone[j + 1L] -
            successful$K_monotone[j]
        )
        fractions <- if (order_gap >= 2L)
          c(0.25, 0.50, 0.75) else 0.50
        if (lo > 0 && hi > 0) {
          new_lambda <- c(
            new_lambda,
            exp((1 - fractions) * log(lo) + fractions * log(hi))
          )
        } else {
          new_lambda <- c(new_lambda, lo + fractions * (hi - lo))
        }
      }
    }
    attained <- successful$K_monotone[
      is.finite(successful$K_monotone)
    ]
    if (nrow(compressed) < 3L ||
        (length(attained) && min(attained) > 2L)) {
      upper <- max(path$lambda)
      if (!is.finite(upper) || upper <= 0) upper <- 0.01
      new_lambda <- c(new_lambda, min(1.75 * upper, 2))
    }
    new_lambda <- unique(new_lambda[is.finite(new_lambda)])
    new_lambda <- new_lambda[
      !vapply(new_lambda, function(z)
        !is.null(records[[key_of(z)]]), logical(1))
    ]
    if (!length(new_lambda)) break
    expansion <- expansion + 1L
    for (lambda in new_lambda)
      fit_missing_lambda(lambda, adaptive_nstart)
    path <- build_path()
    compressed <- compress_h_path(path)
  }

  list(
    success = TRUE, reason = "path generated",
    path = path, compressed = compressed,
    records = records, kde = kde, grid_expansions = expansion
  )
}

# ---------------------------------------------------------------------------
# Likelihood baselines
# ---------------------------------------------------------------------------

loglik_aepd_mixture <- function(x, pars) {
  log_terms <- vapply(seq_along(pars$pi), function(k)
    log(pmax(pars$pi[k], 1e-300)) +
      log_daepd(x, pars$mu[k], pars$sigma[k],
                pars$alpha[k], pars$tau[k]),
    numeric(length(x)))
  row_max <- apply(log_terms, 1L, max)
  sum(row_max + log(rowSums(exp(log_terms - row_max))))
}

optimize_weighted_aepd_loglik <- function(x, weights, start, bounds,
                                          maxeval = 250L) {
  objective <- function(par) {
    mu <- par[1L]
    sigma <- exp(par[2L])
    alpha <- exp(par[3L])
    tau <- stats::plogis(par[4L])
    -sum(weights * log_daepd(x, mu, sigma, alpha, tau))
  }
  par0 <- c(
    start["mu"], log(start["sigma"]), log(start["alpha"]),
    stats::qlogis(start["tau"])
  )
  lower <- c(
    bounds$mu[1L], bounds$log_sigma[1L],
    bounds$log_alpha[1L], bounds$logit_tau[1L]
  )
  upper <- c(
    bounds$mu[2L], bounds$log_sigma[2L],
    bounds$log_alpha[2L], bounds$logit_tau[2L]
  )
  par0 <- pmin(pmax(par0, lower + 1e-10), upper - 1e-10)
  result <- NULL
  if (requireNamespace("nloptr", quietly = TRUE)) {
    result <- tryCatch(nloptr::nloptr(
      par0, objective, lb = lower, ub = upper,
      opts = list(
        algorithm = "NLOPT_LN_BOBYQA",
        xtol_rel = 1e-7, maxeval = maxeval
      )
    ), error = function(e) NULL)
  }
  if (!is.null(result)) {
    solution <- result$solution
  } else {
    result <- tryCatch(stats::optim(
      par0, objective, method = "L-BFGS-B", lower = lower, upper = upper,
      control = list(maxit = maxeval)
    ), error = function(e) NULL)
    solution <- if (is.null(result)) par0 else result$par
  }
  c(
    mu = solution[1L], sigma = exp(solution[2L]),
    alpha = exp(solution[3L]), tau = stats::plogis(solution[4L])
  )
}

em_aepd_mixture <- function(x, K, init, alpha_lower = 1,
                            alpha_upper = 3, tau_lower = 0.01,
                            tau_upper = 0.99, max_iter = 300L,
                            tolerance = 1e-6, maxeval = 250L) {
  bounds <- make_parameter_bounds(
    x, alpha_lower, alpha_upper, tau_lower, tau_upper
  )
  pars <- sanitize_initial(init, K, bounds)
  loglik <- loglik_aepd_mixture(x, pars)
  reason <- "maximum EM iterations reached"

  for (iteration in seq_len(max_iter)) {
    log_terms <- vapply(seq_len(K), function(k)
      log(pmax(pars$pi[k], 1e-300)) +
        log_daepd(x, pars$mu[k], pars$sigma[k],
                  pars$alpha[k], pars$tau[k]),
      numeric(length(x)))
    row_max <- apply(log_terms, 1L, max)
    responsibility <- exp(log_terms - row_max)
    responsibility <- responsibility / rowSums(responsibility)
    nk <- colSums(responsibility)
    candidate <- pars
    candidate$pi <- pmax(nk / length(x), 1e-10)
    candidate$pi <- candidate$pi / sum(candidate$pi)
    for (k in seq_len(K)) {
      if (nk[k] <= 1e-6) next
      start <- c(
        mu = pars$mu[k], sigma = pars$sigma[k],
        alpha = pars$alpha[k], tau = pars$tau[k]
      )
      update <- optimize_weighted_aepd_loglik(
        x, responsibility[, k], start, bounds, maxeval
      )
      old_value <- sum(
        responsibility[, k] *
          log_daepd(x, start["mu"], start["sigma"],
                    start["alpha"], start["tau"])
      )
      new_value <- sum(
        responsibility[, k] *
          log_daepd(x, update["mu"], update["sigma"],
                    update["alpha"], update["tau"])
      )
      if (is.finite(new_value) && new_value + 1e-10 >= old_value) {
        candidate$mu[k] <- update["mu"]
        candidate$sigma[k] <- update["sigma"]
        candidate$alpha[k] <- update["alpha"]
        candidate$tau[k] <- update["tau"]
      }
    }
    candidate <- sort_parameters(candidate)
    new_loglik <- loglik_aepd_mixture(x, candidate)
    if (!is.finite(new_loglik) || new_loglik + 1e-8 < loglik) {
      reason <- "the observed log-likelihood decreased"
      break
    }
    pars <- candidate
    change <- new_loglik - loglik
    loglik <- new_loglik
    if (change <= tolerance * (1 + abs(loglik))) {
      reason <- "EM tolerance satisfied"
      break
    }
  }
  c(pars, list(
    Khat = K, loglik = loglik, df = 5L * K - 1L,
    reason = reason
  ))
}

fit_aepd_mle_hq <- function(x, K_max, nstart = 10L,
                            alpha_lower = 1, alpha_upper = 3,
                            seed = 1L, max_iter = 300L,
                            maxeval = 250L) {
  bounds <- make_parameter_bounds(x, alpha_lower, alpha_upper)
  fits <- vector("list", K_max)
  for (K in seq_len(K_max)) {
    starts <- make_start_list(
      x, K, nstart, bounds, seed = seed + 1000L * K
    )
    candidates <- lapply(starts, function(init)
      tryCatch(em_aepd_mixture(
        x, K, init, alpha_lower, alpha_upper,
        max_iter = max_iter, maxeval = maxeval
      ), error = function(e) NULL))
    good <- which(vapply(candidates, function(z)
      !is.null(z) && is.finite(z$loglik), logical(1)))
    if (length(good)) {
      ll <- vapply(candidates[good], function(z) z$loglik, numeric(1))
      fits[[K]] <- candidates[[good[which.max(ll)]]]
    }
  }
  score <- vapply(fits, function(fit) {
    if (is.null(fit)) return(Inf)
    -2 * fit$loglik + 2 * fit$df * log(log(length(x)))
  }, numeric(1))
  if (!any(is.finite(score))) return(NULL)
  best <- which.min(score)
  selected <- fits[[best]]
  selected$Khat <- best
  selected$HQ <- score[best]
  selected$order_path <- data.frame(
    K = seq_len(K_max),
    loglik = vapply(fits, function(f)
      if (is.null(f)) NA_real_ else f$loglik, numeric(1)),
    df = 5L * seq_len(K_max) - 1L,
    HQ = score
  )
  selected
}

gm_density_grid <- function(grid, fit) {
  pro <- fit$parameters$pro %||% 1
  variance <- fit$parameters$variance$sigmasq %||%
    fit$parameters$variance$sigma2
  sd <- sqrt(if (length(variance) == 1L)
    rep(variance, fit$G) else variance)
  vapply(seq_len(fit$G), function(k)
    pro[k] * stats::dnorm(grid, fit$parameters$mean[k], sd[k]),
    numeric(length(grid))) |>
    rowSums()
}

# ---------------------------------------------------------------------------
# Label alignment and error summaries
# ---------------------------------------------------------------------------

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    rest <- all_permutations(x[-i])
    cbind(x[i], rest)
  }))
}

location_alignment <- function(estimated, truth) {
  K <- length(truth)
  permutations <- all_permutations(seq_len(K))
  loss <- apply(permutations, 1L, function(p)
    sum((estimated[p] - truth)^2))
  permutations[which.min(loss), ]
}

parameter_mse <- function(fit, truth) {
  K <- length(truth$pi)
  if (fit$Khat != K || length(fit$pi) != K)
    return(c(mu = NA, sigma = NA, alpha = NA, tau = NA))
  p <- location_alignment(fit$mu, truth$mu)
  c(
    mu = mean((fit$mu[p] - truth$mu)^2),
    sigma = mean((fit$sigma[p] - truth$sigma)^2),
    alpha = mean((fit$alpha[p] - truth$alpha)^2),
    tau = mean((fit$tau[p] - truth$tau)^2)
  )
}

gm_parameter_mse <- function(fit, truth) {
  K <- length(truth$pi)
  if (is.null(fit) || fit$G != K)
    return(c(mu = NA, sigma = NA, alpha = NA, tau = NA))
  variance <- fit$parameters$variance$sigmasq %||%
    fit$parameters$variance$sigma2
  estimated_sd <- sqrt(if (length(variance) == 1L)
    rep(variance, K) else variance)
  moments <- aepd_component_moments(
    truth$mu, truth$sigma, truth$alpha, truth$tau
  )
  p <- location_alignment(fit$parameters$mean, moments$mean)
  c(
    mu = mean((fit$parameters$mean[p] - moments$mean)^2),
    sigma = mean((estimated_sd[p] - moments$sd)^2),
    alpha = NA, tau = NA
  )
}

truth_evaluation_grid <- function(truth, tail_probability = 1e-9,
                                  grid_n = 4096L) {
  radial <- vapply(truth$alpha, function(alpha)
    stats::qgamma(1 - tail_probability, shape = 1 / alpha)^(1 / alpha),
    numeric(1))
  lower <- min(truth$mu - truth$sigma * radial / (1 - truth$tau))
  upper <- max(truth$mu + truth$sigma * radial / truth$tau)
  grid <- seq(lower, upper, length.out = grid_n)
  w <- trapz_weights(grid)
  density <- mixture_aepd_grid(
    grid, truth$pi, truth$mu, truth$sigma, truth$alpha, truth$tau
  )
  list(
    grid = grid, w = w,
    density = normalize_grid_density(density, w)
  )
}
