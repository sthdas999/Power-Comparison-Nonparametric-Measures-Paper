# ============================================================
# Bootstrap and Jackknife Resampling Study
# Quadratic dependence model:
#       Y = X^2 + epsilon
#       X ~ N(0,1), epsilon ~ N(0, 0.5^2)
#
# Statistics:
#   1. Bergsma-Dassios tau*
#   2. Distance covariance (dCov)
#   3. Graph-based independence statistic
#
# Sample sizes: n = 50, 100, 200
# Bootstrap replications: B = 2000
# Nominal level: alpha = 0.05
# ============================================================

rm(list = ls())

set.seed(2026)

# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------

required_packages <- c(
  "energy",
  "TauStar",
  "igraph"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(energy)
library(TauStar)
library(igraph)

# ------------------------------------------------------------
# Simulation settings
# ------------------------------------------------------------

n.grid <- c(50, 100, 200)

B <- 2000
alpha <- 0.05

# ------------------------------------------------------------
# Quadratic dependence model
# ------------------------------------------------------------

generate_data <- function(n) {
  
  X <- rnorm(n, mean = 0, sd = 1)
  epsilon <- rnorm(n, mean = 0, sd = 0.5)
  
  Y <- X^2 + epsilon
  
  data.frame(
    X = X,
    Y = Y
  )
}

# ------------------------------------------------------------
# Bergsma-Dassios tau* statistic
# ------------------------------------------------------------

tau_star_stat <- function(X, Y) {
  
  out <- tryCatch(
    {
      TauStar::tauStar(X, Y)
    },
    error = function(e) NA_real_
  )
  
  if (length(out) > 1) {
    out <- as.numeric(out[1])
  }
  
  as.numeric(out)
}

# ------------------------------------------------------------
# Distance covariance statistic
# ------------------------------------------------------------

dcov_stat <- function(X, Y) {
  
  out <- tryCatch(
    {
      energy::dcov(X, Y)
    },
    error = function(e) NA_real_
  )
  
  as.numeric(out)
}

# ------------------------------------------------------------
# Graph-based dependence statistic
#
# Here the statistic is constructed from a symmetric
# nearest-neighbour graph. The value is the normalized
# number of graph edges connecting observations that are
# close in the joint (X,Y) space relative to the corresponding
# marginal-neighbour structure.
#
# This implementation is kept explicit so that the graph
# statistic is reproducible without relying on a package-
# specific independence-test wrapper.
# ------------------------------------------------------------

graph_stat <- function(X, Y, k = 1) {
  
  n <- length(X)
  
  if (n < 3) {
    return(NA_real_)
  }
  
  XY <- cbind(
    as.numeric(scale(X)),
    as.numeric(scale(Y))
  )
  
  D <- as.matrix(dist(XY))
  diag(D) <- Inf
  
  # k-nearest-neighbour graph in the joint space
  nn <- apply(
    D,
    1,
    function(z) order(z)[seq_len(min(k, length(z)))]
  )
  
  # Convert nearest-neighbour relationships into an
  # undirected edge set
  edges <- list()
  
  for (i in seq_len(n)) {
    for (j in nn[, i]) {
      
      a <- min(i, j)
      b <- max(i, j)
      
      edges[[length(edges) + 1L]] <- c(a, b)
    }
  }
  
  edges <- unique(do.call(
    rbind,
    edges
  ))
  
  if (nrow(edges) == 0) {
    return(NA_real_)
  }
  
  # Marginal rank-neighbour criterion
  rx <- rank(X, ties.method = "average")
  ry <- rank(Y, ties.method = "average")
  
  concordant_edges <- apply(
    edges,
    1,
    function(e) {
      
      i <- e[1]
      j <- e[2]
      
      dx <- abs(rx[i] - rx[j])
      dy <- abs(ry[i] - ry[j])
      
      # Small rank distance in both margins
      dx <= sqrt(n) && dy <= sqrt(n)
    }
  )
  
  mean(concordant_edges)
}

# ------------------------------------------------------------
# Function to calculate all three statistics
# ------------------------------------------------------------

calculate_statistics <- function(dat) {
  
  X <- dat$X
  Y <- dat$Y
  
  c(
    tau_star = tau_star_stat(X, Y),
    dCov = dcov_stat(X, Y),
    graph = graph_stat(X, Y)
  )
}

# ------------------------------------------------------------
# Bootstrap procedure
# ------------------------------------------------------------

bootstrap_statistic <- function(dat, statistic_function, B = 2000) {
  
  n <- nrow(dat)
  
  boot_values <- numeric(B)
  
  for (b in seq_len(B)) {
    
    ind <- sample.int(
      n = n,
      size = n,
      replace = TRUE
    )
    
    boot_data <- dat[ind, , drop = FALSE]
    
    boot_values[b] <- statistic_function(
      boot_data$X,
      boot_data$Y
    )
  }
  
  boot_values <- boot_values[
    is.finite(boot_values)
  ]
  
  if (length(boot_values) < 2) {
    return(
      list(
        values = boot_values,
        se = NA_real_,
        ci = c(NA_real_, NA_real_)
      )
    )
  }
  
  list(
    values = boot_values,
    se = sd(boot_values),
    ci = as.numeric(
      quantile(
        boot_values,
        probs = c(alpha / 2, 1 - alpha / 2),
        na.rm = TRUE,
        names = FALSE
      )
    )
  )
}

# ------------------------------------------------------------
# Jackknife procedure
# ------------------------------------------------------------

jackknife_statistic <- function(dat, statistic_function) {
  
  n <- nrow(dat)
  
  full_stat <- statistic_function(
    dat$X,
    dat$Y
  )
  
  jack_values <- numeric(n)
  
  for (i in seq_len(n)) {
    
    jack_data <- dat[-i, , drop = FALSE]
    
    jack_values[i] <- statistic_function(
      jack_data$X,
      jack_data$Y
    )
  }
  
  valid <- is.finite(jack_values)
  
  jack_values <- jack_values[valid]
  
  if (length(jack_values) < 2) {
    return(
      list(
        statistic = full_stat,
        values = jack_values,
        pseudovalues = NA_real_,
        se = NA_real_
      )
    )
  }
  
  # Use the number of valid leave-one-out samples
  # for numerical robustness.
  m <- length(jack_values)
  
  pseudovalues <- m * full_stat -
    (m - 1) * jack_values
  
  jack_mean <- mean(pseudovalues)
  
  jack_var <- sum(
    (pseudovalues - jack_mean)^2
  ) / (m * (m - 1))
  
  list(
    statistic = full_stat,
    values = jack_values,
    pseudovalues = pseudovalues,
    se = sqrt(jack_var)
  )
}

# ------------------------------------------------------------
# Wrapper for one statistic
# ------------------------------------------------------------

analyse_one_statistic <- function(
    dat,
    statistic_name,
    statistic_function,
    B = 2000) {
  
  estimate <- statistic_function(
    dat$X,
    dat$Y
  )
  
  boot <- bootstrap_statistic(
    dat = dat,
    statistic_function = statistic_function,
    B = B
  )
  
  jack <- jackknife_statistic(
    dat = dat,
    statistic_function = statistic_function
  )
  
  data.frame(
    Statistic = statistic_name,
    Estimate = estimate,
    Bootstrap_SE = boot$se,
    Jackknife_SE = jack$se,
    CI_Lower = boot$ci[1],
    CI_Upper = boot$ci[2],
    Bootstrap_Replications = length(boot$values),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Run the complete simulation
# ------------------------------------------------------------

results_list <- list()

counter <- 1

for (n in n.grid) {
  
  cat("\n============================================\n")
  cat("Sample size:", n, "\n")
  cat("============================================\n")
  
  dat <- generate_data(n)
  
  # ----------------------------------------------------------
  # tau*
  # ----------------------------------------------------------
  
  cat("Computing tau* ...\n")
  
  results_list[[counter]] <- cbind(
    n = n,
    analyse_one_statistic(
      dat = dat,
      statistic_name = "tau*",
      statistic_function = tau_star_stat,
      B = B
    )
  )
  
  counter <- counter + 1
  
  # ----------------------------------------------------------
  # dCov
  # ----------------------------------------------------------
  
  cat("Computing dCov ...\n")
  
  results_list[[counter]] <- cbind(
    n = n,
    analyse_one_statistic(
      dat = dat,
      statistic_name = "dCov",
      statistic_function = dcov_stat,
      B = B
    )
  )
  
  counter <- counter + 1
  
  # ----------------------------------------------------------
  # Graph-based statistic
  # ----------------------------------------------------------
  
  cat("Computing graph statistic ...\n")
  
  results_list[[counter]] <- cbind(
    n = n,
    analyse_one_statistic(
      dat = dat,
      statistic_name = "Graph-based",
      statistic_function = graph_stat,
      B = B
    )
  )
  
  counter <- counter + 1
}

# ------------------------------------------------------------
# Combine results
# ------------------------------------------------------------

resampling_results <- do.call(
  rbind,
  results_list
)

rownames(resampling_results) <- NULL

# ------------------------------------------------------------
# Print complete results
# ------------------------------------------------------------

print(
  resampling_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# Compact results table
# ------------------------------------------------------------

compact_results <- resampling_results[
  ,
  c(
    "n",
    "Statistic",
    "Estimate",
    "Bootstrap_SE",
    "Jackknife_SE",
    "CI_Lower",
    "CI_Upper"
  )
]

print(
  compact_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# Save numerical results
# ------------------------------------------------------------

write.csv(
  resampling_results,
  file = "Bootstrap_Jackknife_Resampling_Results.csv",
  row.names = FALSE
)

write.csv(
  compact_results,
  file = "Bootstrap_Jackknife_Compact_Results.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Bootstrap vs Jackknife SE comparison
# ------------------------------------------------------------

compact_results$SE_Ratio <-
  compact_results$Jackknife_SE /
  compact_results$Bootstrap_SE

print(
  compact_results[
    ,
    c(
      "n",
      "Statistic",
      "Bootstrap_SE",
      "Jackknife_SE",
      "SE_Ratio"
    )
  ],
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# Optional: visualize bootstrap distributions
# ------------------------------------------------------------

# The following block uses a fresh dataset for each n and
# stores the bootstrap distributions for graphical inspection.

bootstrap_distributions <- list()

for (n in n.grid) {
  
  cat("\nGenerating bootstrap distributions for n =", n, "\n")
  
  dat <- generate_data(n)
  
  boot_tau <- bootstrap_statistic(
    dat,
    tau_star_stat,
    B = B
  )
  
  boot_dcov <- bootstrap_statistic(
    dat,
    dcov_stat,
    B = B
  )
  
  boot_graph <- bootstrap_statistic(
    dat,
    graph_stat,
    B = B
  )
  
  bootstrap_distributions[[paste0("n", n)]] <- list(
    tau_star = boot_tau$values,
    dCov = boot_dcov$values,
    graph = boot_graph$values
  )
}

# ------------------------------------------------------------
# Optional density plots
# ------------------------------------------------------------

for (n in n.grid) {
  
  obj <- bootstrap_distributions[[paste0("n", n)]]
  
  par(mfrow = c(1, 3))
  
  if (length(obj$tau_star) > 1) {
    plot(
      density(obj$tau_star, na.rm = TRUE),
      main = paste("Bootstrap:", expression(tau^"*"),
                   "\nn =", n),
      xlab = expression(tau^"*")
    )
  }
  
  if (length(obj$dCov) > 1) {
    plot(
      density(obj$dCov, na.rm = TRUE),
      main = paste("Bootstrap: dCov\nn =", n),
      xlab = "dCov"
    )
  }
  
  if (length(obj$graph) > 1) {
    plot(
      density(obj$graph, na.rm = TRUE),
      main = paste(
        "Bootstrap: Graph-based\nn =", n
      ),
      xlab = "Graph statistic"
    )
  }
  
  par(mfrow = c(1, 1))
}

# ============================================================
# End of Bootstrap and Jackknife Resampling Study
# ============================================================
