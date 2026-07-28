# Minimal numerical library for the Scenario 3 experiment.

trapz_weights <- function(x) {
  n <- length(x)
  if (n < 2L || any(!is.finite(x)) || any(diff(x) <= 0))
    stop("Quadrature grid must contain at least two increasing values.")
  dx <- diff(x)
  w <- numeric(n)
  w[1L] <- dx[1L] / 2
  w[n] <- dx[n - 1L] / 2
  if (n > 2L)
    w[2L:(n - 1L)] <- (dx[1L:(n - 2L)] + dx[2L:(n - 1L)]) / 2
  w
}

normalize_grid_density <- function(y, w, tolerance = 1e-15) {
  y <- pmax(as.numeric(y), 0)
  mass <- sum(y * w)
  if (!is.finite(mass) || mass <= tolerance)
    stop("Grid density has nonpositive or nonfinite mass.")
  y / mass
}

hellinger2_grid <- function(p, q, w) {
  p <- pmax(p, 0)
  q <- pmax(q, 0)
  2 * pmax(1 - sum(sqrt(p * q) * w), 0)
}

mcp_value <- function(pi, lambda, gamma) {
  sum(ifelse(
    pi <= gamma * lambda,
    lambda * pi - pi^2 / (2 * gamma),
    gamma * lambda^2 / 2
  ))
}

mcp_prime <- function(pi, lambda, gamma)
  pmax(lambda - pi / gamma, 0)

daepd <- function(x, mu, sigma, alpha, tau) {
  sigma <- pmax(sigma, 1e-10)
  alpha <- pmax(alpha, 1e-10)
  tau <- pmin(pmax(tau, 1e-8), 1 - 1e-8)
  constant <- alpha * tau * (1 - tau) /
    (gamma(1 / alpha) * sigma)
  z <- abs(x - mu)^alpha / sigma^alpha
  density <- constant * exp(
    -z * ifelse(x < mu, (1 - tau)^alpha, tau^alpha)
  )
  pmax(density, 1e-300)
}

raepd <- function(n, mu, sigma, alpha, tau) {
  side <- ifelse(stats::runif(n) < (1 - tau), 1, -1)
  radius <- stats::rgamma(n, shape = 1 / alpha)^(1 / alpha)
  side_scale <- ifelse(side == 1, tau, 1 - tau)
  mu + side * sigma * radius / side_scale
}

rmix_aepd <- function(n, pi, mu, sigma, alpha, tau) {
  K <- length(pi)
  counts <- as.vector(stats::rmultinom(1L, n, pi))
  x <- numeric(n)
  position <- 0L
  for (k in seq_len(K)) {
    if (!counts[k]) next
    indices <- position + seq_len(counts[k])
    x[indices] <- raepd(
      counts[k], mu[k], sigma[k], alpha[k], tau[k]
    )
    position <- position + counts[k]
  }
  sample(x)
}

init_params_kmeans <- function(x, K, alpha_init = 2) {
  fit <- stats::kmeans(matrix(x, ncol = 1L), centers = K, nstart = 20L)
  ordering <- order(fit$centers[, 1L])
  cluster <- match(fit$cluster, ordering)
  pi <- tabulate(cluster, nbins = K) / length(x)
  global_sd <- stats::sd(x)
  sigma <- vapply(seq_len(K), function(k) {
    value <- suppressWarnings(stats::sd(x[cluster == k]))
    if (!is.finite(value) || value < 0.1 * global_sd)
      value <- global_sd / 3
    value
  }, numeric(1))
  list(
    pi = pi / sum(pi),
    mu = sort(fit$centers[, 1L]),
    sigma = pmax(sigma, 0.1 * global_sd),
    alpha = rep(alpha_init, K),
    tau = rep(0.5, K)
  )
}

pmhd_fixed_order_fit_one <- function(
    x, K, grid, w, phat, init = NULL,
    tolerance = 1e-5, max_iter = 100L, maxeval = 250L,
    alpha_lower = 1, alpha_upper = 3,
    tau_lower = 0.01, tau_upper = 0.99) {
  sigma_lower <- 1e-3 * stats::sd(x)
  sigma_upper <- 3 * diff(range(x))
  if (is.null(init)) init <- init_params_kmeans(x, K)
  pi <- pmax(init$pi, 1e-8)
  pi <- pi / sum(pi)
  mu <- init$mu
  sigma <- init$sigma
  alpha <- init$alpha
  tau <- init$tau
  previous <- Inf

  for (iteration in seq_len(max_iter)) {
    old_components <- vapply(seq_len(K), function(k)
      daepd(grid, mu[k], sigma[k], alpha[k], tau[k]),
      numeric(length(grid)))
    old_mixture <- pmax(as.vector(old_components %*% pi), 1e-300)
    responsibility <- sweep(
      sweep(old_components, 2L, pi, "*"), 1L, old_mixture, "/"
    )
    normalizer <- colSums(phat * responsibility * w)

    for (k in seq_len(K)) {
      target <- phat * responsibility[, k] / pmax(normalizer[k], 1e-12)
      objective <- function(parameter) {
        candidate <- daepd(
          grid, parameter[1L], exp(parameter[2L]),
          exp(parameter[3L]), stats::plogis(parameter[4L])
        )
        1 - sum(sqrt(target * candidate) * w)
      }
      optimized <- tryCatch(
        nloptr::nloptr(
          x0 = c(mu[k], log(sigma[k]), log(alpha[k]), stats::qlogis(tau[k])),
          eval_f = objective,
          lb = c(
            -Inf, log(sigma_lower), log(alpha_lower),
            stats::qlogis(tau_lower)
          ),
          ub = c(
            Inf, log(sigma_upper), log(alpha_upper),
            stats::qlogis(tau_upper)
          ),
          opts = list(
            algorithm = "NLOPT_LN_BOBYQA",
            xtol_rel = 1e-6, maxeval = maxeval
          )
        ),
        error = function(e) NULL
      )
      if (!is.null(optimized)) {
        mu[k] <- optimized$solution[1L]
        sigma[k] <- exp(optimized$solution[2L])
        alpha[k] <- exp(optimized$solution[3L])
        tau[k] <- stats::plogis(optimized$solution[4L])
      }
    }

    new_components <- vapply(seq_len(K), function(k)
      daepd(grid, mu[k], sigma[k], alpha[k], tau[k]),
      numeric(length(grid)))
    affinity <- vapply(seq_len(K), function(k)
      2 * sum(sqrt(pmax(
        phat * responsibility[, k] * new_components[, k], 0
      )) * w),
      numeric(1))
    # At lambda=0, the simplex update has a closed form.
    pi <- affinity^2
    if (!all(is.finite(pi)) || sum(pi) <= 0) pi <- rep(1 / K, K)
    else pi <- pi / sum(pi)

    mixture <- normalize_grid_density(
      pmax(as.vector(new_components %*% pi), 1e-300), w
    )
    loss <- hellinger2_grid(phat, mixture, w)
    if (abs(previous - loss) < tolerance) break
    previous <- loss
  }

  list(
    pi = pi, mu = mu, sigma = sigma, alpha = alpha, tau = tau,
    Khat = K, H2 = loss, iterations = iteration
  )
}

parameter_names <- function(K) {
  c(
    paste0("pi", seq_len(K - 1L)),
    paste0("mu", seq_len(K)),
    paste0("sigma", seq_len(K)),
    paste0("alpha", seq_len(K)),
    paste0("tau", seq_len(K))
  )
}

pack_theta <- function(parameters) {
  K <- length(parameters$pi)
  value <- c(
    parameters$pi[seq_len(K - 1L)], parameters$mu,
    parameters$sigma, parameters$alpha, parameters$tau
  )
  names(value) <- parameter_names(K)
  value
}

align_by_location <- function(fit) {
  ordering <- order(fit$mu)
  list(
    pi = fit$pi[ordering], mu = fit$mu[ordering],
    sigma = fit$sigma[ordering], alpha = fit$alpha[ordering],
    tau = fit$tau[ordering]
  )
}

next_power_of_two <- function(x)
  2^ceiling(log2(max(x, 2)))

make_theory_kde <- function(
    x, bandwidth, grid_n_min = 2048L, grid_n_max = 262144L,
    points_per_bandwidth = 10) {
  support <- sqrt(5) * bandwidth
  lower <- min(x) - 1.05 * support
  upper <- max(x) + 1.05 * support
  required <- ceiling(
    (upper - lower) / (bandwidth / points_per_bandwidth)
  ) + 1L
  grid_n <- as.integer(max(grid_n_min, next_power_of_two(required)))
  if (grid_n > grid_n_max) {
    stop(
      "KDE requires ", grid_n, " grid points, above limit ", grid_n_max
    )
  }
  kde <- stats::density(
    x, bw = bandwidth, kernel = "epanechnikov",
    n = grid_n, from = lower, to = upper, cut = 0
  )
  w <- trapz_weights(kde$x)
  resolution <- max(diff(kde$x)) / bandwidth
  if (!is.finite(resolution) ||
      resolution > 1 / points_per_bandwidth * 1.001)
    stop("KDE resolution check failed.")
  list(
    grid = kde$x, w = w,
    phat = normalize_grid_density(kde$y, w),
    grid_n = grid_n, resolution = resolution
  )
}

aepd_component_scores <- function(x, mu, sigma, alpha, tau) {
  u <- x - mu
  right <- u >= 0
  q <- ifelse(right, tau, 1 - tau)
  absolute <- abs(u)
  base <- q * absolute / sigma
  z <- base^alpha
  z_log_base <- numeric(length(base))
  positive <- base > 0
  z_log_base[positive] <- z[positive] * log(base[positive])
  cbind(
    mu = alpha * sign(u) * q^alpha * absolute^(alpha - 1) / sigma^alpha,
    sigma = (-1 + alpha * z) / sigma,
    alpha = 1 / alpha + digamma(1 / alpha) / alpha^2 - z_log_base,
    tau = 1 / tau - 1 / (1 - tau) +
      ifelse(right, -alpha * z / tau, alpha * z / (1 - tau))
  )
}

information_matrix <- function(
    truth, tail_probability = 1e-10, grid_n = 200001L) {
  K <- length(truth$pi)
  radial <- vapply(seq_len(K), function(k)
    stats::qgamma(
      1 - tail_probability, shape = 1 / truth$alpha[k]
    )^(1 / truth$alpha[k]),
    numeric(1))
  lower <- min(truth$mu - truth$sigma * radial / (1 - truth$tau))
  upper <- max(truth$mu + truth$sigma * radial / truth$tau)
  x <- sort(unique(c(seq(lower, upper, length.out = grid_n), truth$mu)))
  w <- trapz_weights(x)

  components <- vapply(seq_len(K), function(k)
    daepd(x, truth$mu[k], truth$sigma[k], truth$alpha[k], truth$tau[k]),
    numeric(length(x)))
  weighted <- sweep(components, 2L, truth$pi, "*")
  density <- pmax(rowSums(weighted), 1e-300)
  responsibility <- weighted / density
  score <- matrix(0, nrow = length(x), ncol = 5L * K - 1L)
  colnames(score) <- parameter_names(K)
  score[, seq_len(K - 1L)] <-
    (components[, seq_len(K - 1L), drop = FALSE] - components[, K]) /
    density
  offset <- K - 1L
  for (k in seq_len(K)) {
    component_score <- aepd_component_scores(
      x, truth$mu[k], truth$sigma[k], truth$alpha[k], truth$tau[k]
    )
    score[, offset + k] <- responsibility[, k] * component_score[, "mu"]
    score[, offset + K + k] <-
      responsibility[, k] * component_score[, "sigma"]
    score[, offset + 2L * K + k] <-
      responsibility[, k] * component_score[, "alpha"]
    score[, offset + 3L * K + k] <-
      responsibility[, k] * component_score[, "tau"]
  }
  weighted_score <- score * sqrt(density * w)
  information <- crossprod(weighted_score)
  information <- (information + t(information)) / 2
  attr(information, "captured_mass") <- sum(density * w)
  attr(information, "limits") <- c(lower, upper)
  information
}

safe_solve <- function(matrix, tolerance = 1e-10) {
  decomposition <- eigen((matrix + t(matrix)) / 2, symmetric = TRUE)
  cutoff <- tolerance * max(decomposition$values)
  if (min(decomposition$values) <= cutoff)
    stop("Information matrix is numerically singular.")
  decomposition$vectors %*%
    diag(1 / decomposition$values, nrow = length(decomposition$values)) %*%
    t(decomposition$vectors)
}
