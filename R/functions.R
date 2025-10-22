# =======================================================
# SamplEffCalc R functions - all in one file
# =======================================================
#' SamplEffCalc: Sample Size and Distribution Analysis
#'
#' The SamplEffCalc package provides tools to compute required sample sizes
#' for various distributions, evaluate absolute or relative precision,
#' and visualize observed data with fitted distributions.
#'
#' @docType package
#' @keywords internal
"_PACKAGE"
#'

# -------------------------------------------------------
# Function 1: Histogram with Fitted Distributions (Discrete)
# -------------------------------------------------------
#' Histogram with Fitted Distributions (Discrete)
#'
#' Plots a histogram of observed count data and overlays
#' Poisson, Negative Binomial, Zero-Inflated Poisson (ZIP), and Normal distributions.
#'
#' @param y Numeric vector of counts.
#' @return A ggplot object.
#' @examples
#' y <- rpois(100, 5)
#' plot_distributionsD(y)
#' @export
plot_distributionsD <- function(y) {
  library(ggplot2)

  lambda <- mean(y)
  mu <- lambda
  sd_y <- sd(y)

  bineg <- MASS::glm.nb(y ~ 1)
  theta <- bineg$theta
  size <- theta
  prob <- size / (size + mu)

  x_vals <- min(y):max(y)

  pois_probs <- dpois(x_vals, lambda)
  nbin_probs <- dnbinom(x_vals, size, prob)
  nnorm_probs <- dnorm(x_vals, mean = mu, sd = sd_y)

  p0_obs <- mean(y == 0)
  pi_zip <- max(0, p0_obs - dpois(0, lambda))
  zip_probs <- pi_zip * (x_vals == 0) + (1 - pi_zip) * dpois(x_vals, lambda)

  df_lines <- data.frame(
    x = rep(x_vals, 4),
    y = c(pois_probs, nbin_probs, nnorm_probs, zip_probs),
    dist = factor(rep(c("Poisson", "Negative Binomial", "Normal", "ZIP"), each = length(x_vals)))
  )

  df <- data.frame(y = y)

  ggplot2::ggplot(df, aes(x = y)) +
    geom_histogram(aes(y = after_stat(density)),
                   breaks = seq(min(y)-0.5, max(y)+0.5, by = 1),
                   fill = "lightgray", color = "black") +
    geom_line(data = df_lines, aes(x = x, y = y, color = dist), linewidth = 1) +
    geom_point(data = df_lines, aes(x = x, y = y, color = dist), size = 2) +
    labs(title = "Distribution Histogram",
         subtitle = deparse(substitute(y)),
         x = "Counts",
         y = "Density",
         color = "Distribution") +
    theme_minimal()
}

# -------------------------------------------------------
# Function 2: Histogram with Fitted Distributions (Continuous)
# -------------------------------------------------------
#' Histogram with Fitted Distributions (Continuous)
#'
#' Plots a histogram of continuous or positive data and overlays
#' Normal, Student-t,  Gamma (MLE), and Lognormal (MLE) distributions.
#'
#' @param y Numeric vector of positive values.
#' @return A ggplot object.
#' @examples
#' y <- rgamma(100, shape = 2, rate = 0.5)
#' plot_distributionsC(y)
#' @export
plot_distributionsC <- function(y)
{
  library(ggplot2)
  if (!is.numeric(y))
    stop("y must be numeric")

  mu_norm <- mean(y)
  sd_norm <- sd(y)

  # Fit t-distribution (estimate degrees of freedom using MLE)
  fit_t <- tryCatch({
    fitdistrplus::fitdist(y, "t", start = list(df = length(y) - 1), method = "mle")
  }, error = function(e) {
    # Fallback if MLE fails - use simple method of moments
    # For t-distribution: variance = df/(df-2) for df > 2
    sample_var <- var(y)
    df_est <- ifelse(sample_var > 1, 2 + 2/(sample_var - 1), 3)
    list(estimate = c(df = max(2.1, df_est)))  # df must be > 2
  })
  df_t <- fit_t$estimate["df"]

  y_pos <- y[y > 0]
  if (length(y_pos) < length(y)) {
    warning("Gamma and Lognormal require positive values. Non-positive values ignored.")
  }

  fit_gamma <- fitdistrplus::fitdist(y_pos, "gamma", method = "mle")
  shape_gamma <- fit_gamma$estimate["shape"]
  rate_gamma <- fit_gamma$estimate["rate"]

  fit_lognorm <- fitdistrplus::fitdist(y_pos, "lnorm", method = "mle")
  meanlog_lognorm <- fit_lognorm$estimate["meanlog"]
  sdlog_lognorm <- fit_lognorm$estimate["sdlog"]

  x_vals <- seq(min(y), max(y), length.out = 200)

  dens_norm <- dnorm(x_vals, mean = mu_norm, sd = sd_norm)

  # T-distribution density (using location=mean(y), scale=sd(y), with estimated df)
  dens_t <- dt((x_vals - mean(y))/sd(y), df = df_t) / sd(y)

  dens_gamma <- rep(NA, length(x_vals))
  dens_gamma[x_vals > 0] <- dgamma(x_vals[x_vals > 0], shape = shape_gamma,
                                   rate = rate_gamma)

  dens_lognorm <- rep(NA, length(x_vals))
  dens_lognorm[x_vals > 0] <- dlnorm(x_vals[x_vals > 0], meanlog = meanlog_lognorm,
                                     sdlog = sdlog_lognorm)

  df_lines <- data.frame(x = rep(x_vals, 4),
                         y = c(dens_norm, dens_t, dens_gamma, dens_lognorm),
                         dist = factor(rep(c("Normal", "t-distribution (MLE)",
                                             "Gamma (MLE)", "Lognormal (MLE)"),
                                           each = length(x_vals))))

  df <- data.frame(y = y)
  ggplot2::ggplot(df, aes(x = y)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30, fill = "lightgray", color = "black") +
    geom_line(data = df_lines,
              aes(x = x, y = y, color = dist), linewidth = 1) +
    labs(title = "Continuous Data Histogram with Fitted Distributions",
         subtitle = deparse(substitute(y)),
         x = "Values", y = "Density",
         color = "Distribution") +
    theme_minimal()
}

# -------------------------------------------------------
# Function 3: Relative Error Sample Size
# -------------------------------------------------------
#' Sample Size Calculation for Relative Error
#'
#' Computes sample sizes for both Discrete (Normal, Student-t, Poisson, and Negative Binomial) or Continuous (Normal, Lognormal and Gamma) distributions.
#'
#' @param y Numeric vector of counts.
#' @param epsilon_rel Desired relative error (0.1).
#' @param conf Confidence level (default = 0.95)
#' @param N Population size (optional, for finite population correction)
#' @param Discrete Type of distribution (default = TRUE)
#' @return Prints required sample sizes.
#'
#' @examples
#' y <- rpois(100, 5)
#' relative_sample(y, 0.1)
#' @export
relative_sample <- function(y, epsilon_rel, conf = 0.95, Discrete = TRUE, N = NULL) {

  # Basic statistics
  mu <- mean(y, na.rm = TRUE)
  sigma <- sd(y, na.rm = TRUE)
  lambda <- mu  # for Poisson

  # Negative Binomial fit
  modelo_nb <- MASS::glm.nb(y ~ 1)
  theta <- modelo_nb$theta

  # Gamma fit
  fit_gamma <- fitdistrplus::fitdist(y, "gamma")
  alpha_gamma <- fit_gamma$estimate["shape"]
  beta_gamma  <- fit_gamma$estimate["rate"]

  # Lognormal fit (shift if non-positive values)
  y_pos <- y
  if (any(y <= 0)) {
    shift <- abs(min(y)) + 1e-6
    y_pos <- y + shift
  }
  fit_logn <- fitdistrplus::fitdist(y_pos, "lnorm")
  meanlog  <- fit_logn$estimate["meanlog"]
  sdlog    <- fit_logn$estimate["sdlog"]

  # --- Sample size formulas ---
  sample_size_normal_rel <- function(mu, sigma, epsilon_rel, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf)/2)
    n0 <- (z^2 * sigma^2) / ((epsilon_rel * mu)^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_P_rel <- function(lambda, epsilon_rel, N = NULL) {
    n0 <- 1 / (lambda * epsilon_rel^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_NB_rel <- function(mu, theta, epsilon_rel, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf)/2)
    var_Y <- mu + (mu^2 / theta)
    n0 <- (z^2 * var_Y) / ((epsilon_rel * mu)^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_Gamma_rel <- function(alpha, beta, epsilon_rel, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf)/2)
    mu <- alpha / beta
    var_Y <- alpha / (beta^2)
    n0 <- (z^2 * var_Y) / ((epsilon_rel * mu)^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_Lognorm_rel <- function(meanlog, sdlog, epsilon_rel, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf)/2)
    mu <- exp(meanlog + 0.5 * sdlog^2)
    var_Y <- (exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)
    n0 <- (z^2 * var_Y) / ((epsilon_rel * mu)^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_t_rel <- function(x, epsilon_rel, conf = 0.95, N = NULL) {
    n_pilot <- length(na.omit(x))
    if (n_pilot < 2) stop("Pilot sample must have at least 2 non-missing observations.")
    mu <- mean(x, na.rm = TRUE)
    sigma <- sd(x, na.rm = TRUE)
    t_value <- qt(1 - (1 - conf)/2, df = n_pilot - 1)
    n0 <- ((t_value * sigma) / (epsilon_rel * mu))^2
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  results <- list()

  if (Discrete) {
    results$Normal   <- sample_size_normal_rel(mu, sigma, epsilon_rel, conf, N)
    results$Poisson  <- sample_size_P_rel(lambda, epsilon_rel, N)
    results$NegBinom <- sample_size_NB_rel(mu, theta, epsilon_rel, conf, N)
  } else {
    results$Normal    <- sample_size_normal_rel(mu, sigma, epsilon_rel, conf, N)
    results$t         <- sample_size_t_rel(y, epsilon_rel, conf, N)
    results$Gamma     <- sample_size_Gamma_rel(alpha_gamma, beta_gamma, epsilon_rel, conf, N)
    results$Lognormal <- sample_size_Lognorm_rel(meanlog, sdlog, epsilon_rel, conf, N)
  }

  # --- CONSISTENT OUTPUT FORMAT ---
  cat("=== Relative Error Sample Size Results ===\n")
  cat("Relative Error:", epsilon_rel, "\n")
  cat("Confidence Level:", conf, "\n")
  cat("Distribution Type:", ifelse(Discrete, "Discrete", "Continuous"), "\n")
  cat("Population Size:", ifelse(is.null(N), "Infinite", N), "\n")
  cat("----------------------------------------\n")

  for (nm in names(results)) {
    if (is.list(results[[nm]])) {
      # Finite population case
      cat(nm, "Distribution:\n")
      cat("  Infinite Population (n0):", results[[nm]]$n0, "\n")
      cat("  Finite Population (n_adj):", results[[nm]]$n_adj, "\n")
    } else {
      # Infinite population case
      cat(nm, "Distribution:", results[[nm]], "\n")
    }
  }
  cat("========================================\n")

  invisible(results)
}

# -------------------------------------------------------
# Function 4: Relative Sample Size with ZIP
# -------------------------------------------------------
#' Sample Size Calculation for Relative Error with ZIP
#'
#' Computes sample sizes for both Discrete (Normal, Student-t, Poisson, and Negative Binomial) or
#' Continuous (Normal, Lognormal and Gamma) distributions.
#'
#' @param y Numeric vector of counts.
#' @param epsilon_rel Desired relative error (0.1)
#' @param conf Confidence level (default 0.95)
#' @param N Population size (optional, for finite population correction)
#' @return Prints required sample sizes.
#' @examples
#' y <- rpois(100, 3); y[sample(1:100, 20)] <- 0
#' relative_sample_zip(y, 0.1)
#' @export
relative_sample_zip <- function(y, epsilon_rel, conf = 0.95, N = NULL) {
  # Load required packages
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required.")

  modelo_nb <- MASS::glm.nb(y ~ 1)
  theta <- modelo_nb$theta
  mu <- mean(y)
  sigma <- sd(y)
  lambda <- mean(y)

  sample_size_NB_ci_Rel <- function(mu, theta, epsilon_rel, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf) / 2)
    var_Y <- mu + (mu^2 / theta)
    n0 <- (z^2 * var_Y) / (epsilon_rel^2 * mu^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_P_ci_Rel <- function(lambda, epsilon_rel, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf) / 2)
    n0 <- (z^2) / (epsilon_rel^2 * lambda)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_normal_rel <- function(mu, sigma, epsilon_rel, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf) / 2)
    epsilon <- epsilon_rel * mu
    n0 <- (z * sigma / epsilon)^2
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_ZIP_rel <- function(y, epsilon_rel, conf = 0.95, N = NULL) {
    estimate_pi_zip_mom <- function(y, lower_mu = 1e-06, upper_mu = 1e3) {
      p0 <- mean(y == 0); ybar <- mean(y)
      g <- function(mu) {
        if(mu <= 0) return(1e6)
        val <- ((1 - p0 + exp(-mu)) / (1 - exp(-mu))) * mu - ybar
        return(val)
      }
      lo <- lower_mu; hi <- upper_mu; found <- FALSE
      for(i in 1:50){
        fa <- g(lo); fb <- g(hi)
        if(is.finite(fa) && is.finite(fb) && fa*fb < 0){
          found <- TRUE; break
        }
        lo <- lo/2; hi <- hi*2
      }
      if(!found) return(list(mu = ybar, pi = 0, p0 = p0, ybar = ybar))
      sol <- uniroot(g, lower = lo, upper = hi)
      mu_hat <- sol$root
      pi_hat <- (p0 - exp(-mu_hat)) / (1 - exp(-mu_hat))
      return(list(mu = mu_hat, pi = pi_hat, p0 = p0, ybar = ybar))
    }
    res <- estimate_pi_zip_mom(y)
    pi <- res$pi; mu <- res$mu
    mu_m <- (1-pi)*mu
    VarY <- (1-pi)*mu + (mu^2)*pi*(1-pi)
    z <- qnorm(1 - (1 - conf)/2)
    n0 <- (z^2 * VarY) / (epsilon_rel^2 * mu_m^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  n1 <- sample_size_NB_ci_Rel(mu, theta, epsilon_rel, conf, N)
  n2 <- sample_size_P_ci_Rel(lambda, epsilon_rel, conf, N)
  n3 <- sample_size_normal_rel(mu, sigma, epsilon_rel, conf, N)
  n4 <- sample_size_ZIP_rel(y, epsilon_rel, conf, N)


  cat("=== Relative Error Sample Size Results (ZIP) ===\n")
  cat("Relative Error:", epsilon_rel, "\n")
  cat("Confidence Level:", conf, "\n")
  cat("Population Size:", ifelse(is.null(N), "Infinite", N), "\n")
  cat("----------------------------------------\n")

  results <- list(
    Negative_Binomial = n1,
    Poisson = n2,
    Normal = n3,
    ZIP = n4
  )

  for (nm in names(results)) {
    if (is.list(results[[nm]])) {
      # Finite population case
      cat(nm, "Distribution:\n")
      cat("  Infinite Population (n0):", results[[nm]]$n0, "\n")
      cat("  Finite Population (n_adj):", results[[nm]]$n_adj, "\n")
    } else {
      # Infinite population case
      cat(nm, "Distribution:", results[[nm]], "\n")
    }
  }
  cat("========================================\n")

  invisible(results)
}

# -------------------------------------------------------
# Function 5: Absolute Error Sample Size
# -------------------------------------------------------
#' Sample Size Calculation for Absolute Error
#'
#' Computes sample sizes for Normal, Poisson, Negative Binomial,
#' Gamma, and Lognormal distributions.
#'
#' @param y Numeric vector of values.
#' @param epsilon Desired absolute error.
#' @param conf Confidence level (default 0.95).
#' @param Discrete If TRUE, use discrete distributions; if FALSE, use continuous.
#' @param N Population size (optional, for finite population correction)
#' @return Prints required sample sizes.
#' @examples
#' y <- rpois(100, 5)
#' absolute_sample(y, 2)
#' @export
absolute_sample <- function(y, epsilon, conf = 0.95, Discrete = TRUE, N = NULL) {

  mu <- mean(y, na.rm = TRUE)
  sigma <- sd(y, na.rm = TRUE)
  lambda <- mu

  # Negative Binomial fit
  modelo_nb <- MASS::glm.nb(y ~ 1)
  theta <- modelo_nb$theta

  # Gamma fit
  fit_gamma <- fitdistrplus::fitdist(y, "gamma")
  alpha_gamma <- fit_gamma$estimate["shape"]
  beta_gamma  <- fit_gamma$estimate["rate"]

  # Lognormal fit
  y_pos <- y
  if (any(y <= 0)) {
    shift <- abs(min(y)) + 1e-6
    y_pos <- y + shift
  }
  fit_logn <- fitdistrplus::fitdist(y_pos, "lnorm")
  meanlog  <- fit_logn$estimate["meanlog"]
  sdlog    <- fit_logn$estimate["sdlog"]

  # Internal functions with finite population correction
  sample_size_normal_abs <- function(mu, sigma, epsilon, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf) / 2)
    n0 <- (z * sigma / epsilon)^2
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_P_abs <- function(lambda, epsilon, N = NULL) {
    n0 <- lambda / (epsilon^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_NB_abs <- function(mu, theta, epsilon, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf) / 2)
    var_Y <- mu + (mu^2 / theta)
    n0 <- z^2 * var_Y / (epsilon^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_Gamma_abs <- function(alpha, beta, epsilon, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf) / 2)
    var_Y <- alpha / (beta^2)
    n0 <- z^2 * var_Y / (epsilon^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_Lognorm_abs <- function(meanlog, sdlog, epsilon, conf = 0.95, N = NULL) {
    z <- qnorm(1 - (1 - conf) / 2)
    mu <- exp(meanlog + 0.5 * sdlog^2)
    var_Y <- (exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2)
    n0 <- z^2 * var_Y / (epsilon^2)
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  sample_size_t_abs <- function(x, epsilon, conf = 0.95, N = NULL) {
    n_pilot <- length(na.omit(x))
    if (n_pilot < 2) stop("Pilot sample must have at least 2 non-missing observations.")
    mu <- mean(x, na.rm = TRUE)
    sigma <- sd(x, na.rm = TRUE)
    t_value <- qt(1 - (1 - conf)/2, df = n_pilot - 1)
    n0 <- (t_value * sigma / epsilon)^2
    if (!is.null(N)) {
      n_adj <- n0 / (1 + (n0 - 1) / N)
      return(list(n0 = ceiling(n0), n_adj = ceiling(n_adj)))
    } else {
      return(ceiling(n0))
    }
  }

  # Results
  results <- list()

  if (Discrete) {
    results$Normal   <- sample_size_normal_abs(mu, sigma, epsilon, conf, N)
    results$Poisson  <- sample_size_P_abs(lambda, epsilon, N)
    results$NegBinom <- sample_size_NB_abs(mu, theta, epsilon, conf, N)
  } else {
    results$Normal    <- sample_size_normal_abs(mu, sigma, epsilon, conf, N)
    results$Student_t <- sample_size_t_abs(y, epsilon, conf, N)
    results$Gamma     <- sample_size_Gamma_abs(alpha_gamma, beta_gamma, epsilon, conf, N)
    results$Lognormal <- sample_size_Lognorm_abs(meanlog, sdlog, epsilon, conf, N)
  }

  cat("=== Absolute Error Sample Size Results ===\n")
  cat("Absolute Error:", epsilon, "\n")
  cat("Confidence Level:", conf, "\n")
  cat("Distribution Type:", ifelse(Discrete, "Discrete", "Continuous"), "\n")
  cat("Population Size:", ifelse(is.null(N), "Infinite", N), "\n")
  cat("----------------------------------------\n")

  for (nm in names(results)) {
    if (is.list(results[[nm]])) {
      # Finite population case
      cat(nm, "Distribution:\n")
      cat("  Infinite Population (n0):", results[[nm]]$n0, "\n")
      cat("  Finite Population (n_adj):", results[[nm]]$n_adj, "\n")
    } else {
      # Infinite population case
      cat(nm, "Distribution:", results[[nm]], "\n")
    }
  }
  cat("========================================\n")

  invisible(results)
}
# -----------------------------------------------
# Function 6: Bootstrap-based Relative Precision
# -----------------------------------------------
#' Bootstrap-based Relative Precision Analysis
#'
#' Estimates sampling distribution of the mean using bootstrap, compute confidence intervals,
#' and track relative precision across increasing sample sizes. It calculates the precision for both an
#' infinite population size or from a complete sample (finite population size).
#' @param y Numeric vector of observations.
#' @param epsilon_rel Numeric. Desired relative precision threshold (e.g. 0.2).
#' @param epsilon Numeric. Desired absolute precision threshold (e.g 1)
#' @param n_boot Number of bootstrap replicates (default = 1000).
#' @param conf Confidence level (default = 0.95).
#' @param seed Random seed (default = 123).
#' @param N_pop Finite population size (default = NULL).
#' @return List with table of results and first sample size reaching desired precision.
#' @examples
#' y <- rpois(30, 5)
#' bootstrap_precision(y, epsilon_rel = 0.2, n_boot = 200)
#' bootstrap_precision(y, epsilon_rel = 0.2, n_boot = 200, N_pop=100)
#' @export
bootstrap_precision <- function(y, epsilon_rel = NULL, epsilon = NULL,
                                 n_boot = 1000, conf = 0.95, seed = 123,
                                 N_pop = NULL) {

  library(ggplot2)
  library(gridExtra)

  set.seed(seed)
  y <- na.omit(y)
  N <- length(y)
  alpha <- 1 - conf
  results <- data.frame(n = integer(), mean = numeric(), lower = numeric(),
                        upper = numeric(), ci_width = numeric(), precision = numeric())

  for (n in 1:N) {
    boot_means <- replicate(n_boot, mean(sample(y, n, replace = TRUE)))
    mean_hat <- mean(boot_means)
    lower <- quantile(boot_means, probs = alpha/2)
    upper <- quantile(boot_means, probs = 1 - alpha/2)
    ci_width <- upper - lower

    # Compute precision: either relative or absolute
    if (!is.null(epsilon_rel)) {
      precision <- (ci_width / abs(mean_hat)) / 2
    } else if (!is.null(epsilon)) {
      precision <- ci_width / 2
    } else {
      stop("Either epsilon_rel or epsilon must be provided.")
    }

    results <- rbind(results, data.frame(n = n, mean = mean_hat,
                                         lower = lower, upper = upper,
                                         ci_width = ci_width, precision = precision))
  }

  # Select first n meeting the precision criterion
  idx <- which(results$precision < ifelse(!is.null(epsilon_rel), epsilon_rel, epsilon))
  if (length(idx) > 0) {
    n0 <- results$n[idx[1]]
    ci_info <- results[results$n == n0, c("mean", "lower", "upper", "ci_width", "precision")]
  } else {
    warning("Desired precision not achieved with available sample size.")
    return(NA)
  }

  # Finite population adjustment
  if (!is.null(N_pop)) {
    n_adj <- n0 / (1 + (n0 - 1) / N_pop)
    n_adj <- ceiling(n_adj)
    sample_sizes <- list(n0 = n0, n_adj = n_adj)
  } else {
    sample_sizes <- list(n0 = n0, n_adj = NULL)
  }

  # Plotting
  p1 <- ggplot2::ggplot(results, aes(x = n, y = mean)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
    geom_vline(xintercept = n0, linetype = "dashed", color = "green") +
    geom_hline(yintercept = ci_info$mean, linetype = "dashed", color = "gray") +
    labs(title = "Bootstrap Mean and Confidence Intervals", x = "Sample Size (n)", y = "Mean ± CI") +
    annotate("text", x = n0, y = min(results$mean), label = paste("n0 =", n0), vjust = -1, color = "green") +
    theme_minimal()

  p2 <- ggplot2::ggplot(results, aes(x = n, y = precision)) +
    geom_line(color = "red") +
    geom_hline(yintercept = ifelse(!is.null(epsilon_rel), epsilon_rel, epsilon),
               linetype = "dashed", color = "blue") +
    geom_vline(xintercept = n0, linetype = "dashed", color = "green") +
    labs(title = "Precision vs Sample Size", x = "Sample Size (n)", y = ifelse(!is.null(epsilon_rel), "Relative Precision", "Absolute Precision")) +
    annotate("text", x = n0, y = max(results$precision), label = paste("n0 =", n0), vjust = -1, color = "green") +
    theme_minimal()

  if (!is.null(N_pop)) {
    p1 <- p1 + geom_vline(xintercept = n_adj, linetype = "dashed", color = "purple") +
      annotate("text", x = n_adj, y = max(results$mean), label = paste("n_adj =", n_adj), vjust = -1, color = "purple")
    p2 <- p2 + geom_vline(xintercept = n_adj, linetype = "dashed", color = "purple") +
      annotate("text", x = n_adj, y = max(results$precision) * 0.9, label = paste("n_adj =", n_adj), vjust = -1, color = "purple")
  }

  gridExtra::grid.arrange(p1, p2, nrow = 2)


  cat("=== Bootstrap Precision Sample Size Results ===\n")
  cat("Required Sample Size (n0):", n0, "\n")
  if (!is.null(N_pop)) {
    cat("Adjusted Sample Size (n_adj):", n_adj, "\n")
    cat("Population Size (N):", N_pop, "\n")
  } else {
    cat("Population Size: Infinite\n")
  }
  cat("Precision Target:", ifelse(!is.null(epsilon_rel), epsilon_rel, epsilon), "\n")
  cat("Achieved Precision:", round(ci_info$precision, 4), "\n")
  cat("Confidence Level:", conf, "\n")
  cat("Bootstrap Replicates:", n_boot, "\n")
  cat("Estimated Mean:", round(ci_info$mean, 4), "\n")
  cat("Confidence Interval: [", round(ci_info$lower, 4), ", ", round(ci_info$upper, 4), "]\n", sep = "")
  cat("Confidence Interval Width:", round(ci_info$ci_width, 4), "\n")
  cat("Pilot Sample Size:", N, "\n")
  cat("===============================================\n")

  # Return results invisibly
  results_list <- list(
    sample_sizes = sample_sizes,
    infinite_population = list(n0 = n0, mean = ci_info$mean, lower_CI = ci_info$lower,
                               upper_CI = ci_info$upper, ci_width = ci_info$ci_width,
                               precision = ci_info$precision),
    finite_population = if (!is.null(N_pop)) list(n_adj = n_adj, population_size = N_pop) else NULL,
    all_results = results
  )
  invisible(results_list)
}
# -------------------------------------------------------
# Function 7: Finite Population Correction
# -------------------------------------------------------
#' Finite Population Correction for Sample Size
#'
#' Adjusts an initial sample size calculation (n0) for a finite population.
#'
#' @param n0 Numeric. Initial sample size assuming an infinite population.
#' @param N Numeric. Total population size.
#'
#' @return Adjusted sample size (numeric).
#' @examples
#' n0 <- 50     # initial sample size
#' N  <- 200    # total population size
#' finite_population_correction(n0, N)
#' @export
finite_population_correction <- function(n0, N) {
  n_adj <- n0 / (1 + (n0 - 1) / N)
  return(ceiling(n_adj))  # round up to ensure sufficiency
}

# -------------------------------------------------------
# Function 8: Sample Size for Binomial Proportion
# -------------------------------------------------------
#' Sample Size Calculation for a Binomial Proportion
#'
#' Computes the required sample size to estimate a proportion in a binomial
#' distribution with a given absolute or relative error and confidence level.
#'
#' @param x Numeric vector of 0/1 data, or NULL if user provides p.
#' @param p Numeric. Optional expected proportion of "successes" (0 < p < 1). If NULL, estimated from x.
#' @param epsilon Numeric. Desired precision (margin of error, 0 < E < 1).
#' @param conf Numeric. Confidence level (default = 0.95).
#' @param N Optional. Population size for finite population correction (default = NULL).
#' @param relative Logical. If TRUE, interprets E as relative error (default = FALSE).
#'
#' @return Required sample size (integer).
#' @examples
#' # Example with observed data and absolute precision
#' y <- c(rep(1, 30), rep(0, 70))
#' sample_size_binomial(x = y, epsilon = 0.05)
#'
#' # Example with relative precision (5% of the proportion)
#' sample_size_binomial(p = 0.4, epsilon = 0.05, conf = 0.95, N = 500, relative = TRUE)
#' @export
sample_size_binomial <- function(x = NULL, p = NULL, epsilon, conf = 0.95, N = NULL, relative = FALSE) {

  # Validate inputs
  if (missing(epsilon) || epsilon <= 0) stop("Error 'epsilon' must be > 0")
  if (conf <= 0 || conf >= 1) stop("'conf' must be in (0,1)")

  # Estimate p from data if not provided
  if (is.null(p)) {
    if (is.null(x)) stop("Provide either 'x' or 'p'.")
    if (!all(x %in% c(0,1))) stop("x must contain only 0s and 1s.")
    p <- mean(x)
  }

  if (p <= 0 || p >= 1) stop("Proportion 'p' must be between 0 and 1.")

  # Interpret E as absolute or relative
  E_abs <- if (relative) epsilon * p else epsilon

  # Validate absolute error doesn't exceed bounds
  if (E_abs <= 0 || E_abs >= 1) {
    stop("Resulting absolute error must be between 0 and 1")
  }

  # Z-score for confidence level
  z <- qnorm(1 - (1 - conf) / 2)

  # Initial sample size (infinite population)
  n0 <- (z^2 * p * (1 - p)) / (E_abs^2)

  # Finite population correction
  if (!is.null(N)) {
    if (N <= 0) stop("'N' must be > 0")
    n0 <- n0 / (1 + (n0 - 1) / N)
  }

  # Create results list
  results <- list(
    SampleSize = ceiling(n0),
    Proportion = p,
    PrecisionType = ifelse(relative, "Relative", "Absolute"),
    PrecisionUsed = epsilon,
    AbsoluteErrorUsed = E_abs,
    Confidence = conf,
    PopulationSize = ifelse(is.null(N), "Inf", N),
    ConfidenceInterval = c(p - E_abs, p + E_abs)
  )

  # Formatted output
  cat("=== Binomial Proportion Sample Size Results ===\n")
  cat("Sample Size:", results$SampleSize, "\n")
  cat("Proportion:", round(results$Proportion, 4), "\n")
  cat("Precision Type:", results$PrecisionType, "\n")
  cat("Precision Used:", results$PrecisionUsed, "\n")
  cat("Absolute Error:", round(results$AbsoluteErrorUsed, 4), "\n")
  cat("Confidence Level:", results$Confidence, "\n")
  cat("Population Size:", results$PopulationSize, "\n")
  cat("Confidence Interval: [",
      round(results$ConfidenceInterval[1], 4), ", ",
      round(results$ConfidenceInterval[2], 4), "]\n", sep = "")

  invisible(results)
}
# -------------------------------------------------------
# Function 9: Sample Size for Beta Distribution (with 0/1 adjustment)
# -------------------------------------------------------
#' Sample Size Calculation for a Beta Distribution Mean
#'
#' Computes the required sample size to estimate the mean of a Beta distribution
#' with a given absolute precision (margin of error) and confidence level.
#' Automatically transforms any 0/1 values in the data using Smithson & Verkuilen (2006).
#'
#' @param x Numeric vector of data in [0,1], or NULL if user provides alpha and beta.
#' @param alpha Numeric. Shape1 parameter of the Beta distribution (optional if x given).
#' @param beta Numeric. Shape2 parameter of the Beta distribution (optional if x given).
#' @param epsilon Numeric. Desired absolute precision (margin of error, >0).
#' @param conf Numeric. Confidence level (default = 0.95).
#' @param N Optional. Population size for finite population correction (default = NULL).
#' @param relative Logical. If TRUE, interprets E as relative error (default = FALSE).
#'
#' @return Required sample size (integer).
#' @examples
#' # Example with observed data containing 0
#' x <- c(0.2, 0.1, 0, 0, 0, 0.4)
#' sample_size_beta(x = x, epsilon = 0.05)
#'
#' # Example with provided parameters
#' sample_size_beta(alpha = 2, beta = 5, epsilon = 0.05, conf = 0.95, N = 500)
#' @export
sample_size_beta <- function(x = NULL, alpha = NULL, beta = NULL,
                             epsilon, conf = 0.95, N = NULL, relative = FALSE) {

  if (missing(epsilon) || epsilon <= 0) stop("'E' must be > 0")
  if (conf <= 0 || conf >= 1) stop("'conf' must be in (0,1)")

  # Se alpha/beta não fornecidos, estimar a partir de x
  if (is.null(alpha) || is.null(beta)) {
    if (is.null(x)) stop("Provide either 'x' or both 'alpha' and 'beta'.")
    if (!is.numeric(x) || any(x < 0 | x > 1)) stop("'x' must be numeric in [0,1].")

    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 2) stop("At least 2 non-missing observations required.")

    # Smithson & Verkuilen (2006) transformation
    if (any(x == 0 | x == 1)) {
      x <- (x * (n - 1) + 0.5) / n
    }

    # Estimativa por método dos momentos
    m <- mean(x)
    s2 <- var(x)
    if (s2 <= 0 || s2 >= m * (1 - m)) {
      warning("Variance too small or too large; using fallback alpha=1, beta=1")
      alpha <- 1
      beta <- 1
    } else {
      tmp <- (m * (1 - m) / s2) - 1
      alpha <- m * tmp
      beta  <- (1 - m) * tmp
    }
  }

  if (!(alpha > 0 && beta > 0)) stop("alpha and beta must be > 0")


  mu <- alpha / (alpha + beta)
  var_Y <- (alpha * beta) / (((alpha + beta)^2) * (alpha + beta + 1))


  E_abs <- if (relative) epsilon * mu else epsilon

  # Z-score
  z <- qnorm(1 - (1 - conf) / 2)
  n0 <- (z^2 * var_Y) / (E_abs^2)

  # finite population correction
  if (!is.null(N)) {
    if (N <= 0) stop("'N' must be > 0")
    n0 <- n0 / (1 + (n0 - 1) / N)
  }

  # Criar lista de resultados
  results <- list(
    SampleSize = ceiling(n0),
    Mean = mu,
    PrecisionType = ifelse(relative, "Relative", "Absolute"),
    PrecisionUsed = epsilon,
    AbsoluteErrorUsed = E_abs,
    Confidence = conf,
    PopulationSize = ifelse(is.null(N), "Inf", N),
    ConfidenceInterval = c(mu - E_abs, mu + E_abs)
  )

  # Saída formatada (similar à Function 8)
  cat("=== Beta Distribution Sample Size Results ===\n")
  cat("Sample Size:", results$SampleSize, "\n")
  cat("Mean:", round(results$Mean, 4), "\n")
  cat("Precision Type:", results$PrecisionType, "\n")
  cat("Precision Used:", results$PrecisionUsed, "\n")
  cat("Absolute Error:", round(results$AbsoluteErrorUsed, 4), "\n")
  cat("Confidence Level:", results$Confidence, "\n")
  cat("Population Size:", results$PopulationSize, "\n")
  cat("Confidence Interval: [",
      round(results$ConfidenceInterval[1], 4), ", ",
      round(results$ConfidenceInterval[2], 4), "]\n", sep = "")

  invisible(results)
}

