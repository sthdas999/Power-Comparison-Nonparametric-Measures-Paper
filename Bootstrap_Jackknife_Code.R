# ============================================================
# Bootstrap and Jackknife Resampling Study
# Quadratic dependence model:
#       Y = X^2 + epsilon
#       X ~ N(0,1), epsilon ~ N(0, 0.5^2)
#
# Statistics:
#   1. Bergsma-Dassios tau*
#   2. Distance covariance (dCov)
#   3. Graph-based dependence statistic
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
  "TauStar"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(energy)
library(TauStar)

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
  
  X <- rnorm(
    n,
    mean = 0,
    sd = 1
  )
  
  epsilon <- rnorm(
    n,
    mean = 0,
    sd = 0.5
  )
  
  Y <- X^2 + epsilon
  
  data.frame(
    X = X,
    Y = Y
  )
}

# ------------------------------------------------------------
# Bergsma-Dassios tau* statistic
#
# TauStar package uses tStar(), not tauStar().
# ------------------------------------------------------------

tau_star_stat <- function(X, Y) {
  
  out <- tryCatch(
    {
      TauStar::tStar(X, Y)
    },
    error = function(e) {
      NA_real_
    }
  )
  
  as.numeric(out)
}

# ------------------------------------------------------------
# Distance covariance
# ------------------------------------------------------------

dcov_stat <- function(X, Y) {
  
  out <- tryCatch(
    {
      energy::dcov(X, Y)
    },
    error = function(e) {
      NA_real_
    }
  )
  
  as.numeric(out)
}

# ------------------------------------------------------------
# Graph-based dependence statistic
#
# A 1-nearest-neighbour graph is constructed in the
# standardized joint (X,Y) space.
#
# The statistic is the proportion of graph edges for which
# both marginal rank distances are relatively small.
# ------------------------------------------------------------

graph_stat <- function(X, Y, k = 1) {
  
  n <- length(X)
  
  if (n < 3) {
    return(NA_real_)
  }
  
  # Standardized joint coordinates
  XY <- cbind(
    as.numeric(scale(X)),
    as.numeric(scale(Y))
  )
  
  # Pairwise distances
  D <- as.matrix(dist(XY))
  
  diag(D) <- Inf
  
  # ----------------------------------------------------------
  # Obtain k nearest neighbours for every observation.
  #
  # IMPORTANT:
  # Do not use apply() here because k = 1 causes
  # dimension dropping.
  # ----------------------------------------------------------
  
  nn <- lapply(
    seq_len(n),
    function(i) {
      order(D[i, ])[seq_len(k)]
    }
  )
  
  # ----------------------------------------------------------
  # Construct undirected edge list
  # ----------------------------------------------------------
  
  edges <- vector(
    mode = "list",
    length = n * k
  )
  
  counter <- 1L
  
  for (i in seq_len(n)) {
    
    for (j in nn[[i]]) {
      
      a <- min(i, j)
      b <- max(i, j)
      
      edges[[counter]] <- c(a, b)
      
      counter <- counter + 1L
    }
  }
  
  edges <- do.call(
    rbind,
    edges
  )
  
  # Remove duplicate undirected edges
  edges <- unique(edges)
  
  if (nrow(edges) == 0) {
    return(NA_real_)
  }
  
  # ----------------------------------------------------------
  # Marginal ranks
  # ----------------------------------------------------------
  
  rx <- rank(
    X,
    ties.method = "average"
  )
  
  ry <- rank(
    Y,
    ties.method = "average"
  )
  
  # ----------------------------------------------------------
  # Concordance criterion
  # ----------------------------------------------------------
  
  concordant_edges <- apply(
    edges,
    1,
    function(e) {
      
      i <- e[1]
      j <- e[2]
      
      dx <- abs(rx[i] - rx[j])
      dy <- abs(ry[i] - ry[j])
      
      dx <= sqrt(n) &&
        dy <= sqrt(n)
    }
  )
  
  mean(concordant_edges)
}

# ------------------------------------------------------------
# Calculate all statistics
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

bootstrap_statistic <- function(
    dat,
    statistic_function,
    B = 2000) {
  
  n <- nrow(dat)
  
  boot_values <- numeric(B)
  
  for (b in seq_len(B)) {
    
    ind <- sample.int(
      n = n,
      size = n,
      replace = TRUE
    )
    
    boot_data <- dat[
      ind,
      ,
      drop = FALSE
    ]
    
    boot_values[b] <- statistic_function(
      boot_data$X,
      boot_data$Y
    )
  }
  
  # Remove failed evaluations
  boot_values <- boot_values[
    is.finite(boot_values)
  ]
  
  if (length(boot_values) < 2) {
    
    return(
      list(
        values = boot_values,
        se = NA_real_,
        ci = c(
          NA_real_,
          NA_real_
        )
      )
    )
  }
  
  boot_se <- sd(
    boot_values
  )
  
  boot_ci <- quantile(
    boot_values,
    probs = c(
      alpha / 2,
      1 - alpha / 2
    ),
    na.rm = TRUE,
    names = FALSE
  )
  
  list(
    values = boot_values,
    se = boot_se,
    ci = as.numeric(boot_ci)
  )
}

# ------------------------------------------------------------
# Jackknife procedure
# ------------------------------------------------------------

jackknife_statistic <- function(
    dat,
    statistic_function) {
  
  n <- nrow(dat)
  
  # Statistic from complete sample
  full_stat <- statistic_function(
    dat$X,
    dat$Y
  )
  
  # Leave-one-out statistics
  jack_values <- numeric(n)
  
  for (i in seq_len(n)) {
    
    jack_data <- dat[
      -i,
      ,
      drop = FALSE
    ]
    
    jack_values[i] <- statistic_function(
      jack_data$X,
      jack_data$Y
    )
  }
  
  valid <- is.finite(
    jack_values
  )
  
  # If any leave-one-out statistic fails,
  # return NA rather than silently changing n.
  if (!all(valid)) {
    
    return(
      list(
        statistic = full_stat,
        values = jack_values,
        pseudovalues = rep(
          NA_real_,
          n
        ),
        se = NA_real_
      )
    )
  }
  
  # ----------------------------------------------------------
  # Jackknife pseudovalues
  # ----------------------------------------------------------
  
  pseudovalues <-
    n * full_stat -
    (n - 1) * jack_values
  
  jack_mean <- mean(
    pseudovalues
  )
  
  # Jackknife variance
  jack_var <-
    sum(
      (pseudovalues - jack_mean)^2
    ) /
    (n * (n - 1))
  
  list(
    statistic = full_stat,
    values = jack_values,
    pseudovalues = pseudovalues,
    se = sqrt(jack_var)
  )
}

# ------------------------------------------------------------
# Analyse one statistic
# ------------------------------------------------------------

analyse_one_statistic <- function(
    dat,
    statistic_name,
    statistic_function,
    B = 2000) {
  
  # Original estimate
  estimate <- statistic_function(
    dat$X,
    dat$Y
  )
  
  # Bootstrap
  boot <- bootstrap_statistic(
    dat = dat,
    statistic_function = statistic_function,
    B = B
  )
  
  # Jackknife
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
    Bootstrap_Replications =
      length(boot$values),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Main simulation
# ============================================================

results_list <- list()

counter <- 1L

for (n in n.grid) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Sample size:",
    n,
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  # Generate one dataset for this sample size
  dat <- generate_data(n)
  
  # ----------------------------------------------------------
  # tau*
  # ----------------------------------------------------------
  
  cat(
    "Computing tau* ...\n"
  )
  
  results_list[[counter]] <- cbind(
    n = n,
    analyse_one_statistic(
      dat = dat,
      statistic_name = "tau*",
      statistic_function = tau_star_stat,
      B = B
    )
  )
  
  counter <- counter + 1L
  
  # ----------------------------------------------------------
  # dCov
  # ----------------------------------------------------------
  
  cat(
    "Computing dCov ...\n"
  )
  
  results_list[[counter]] <- cbind(
    n = n,
    analyse_one_statistic(
      dat = dat,
      statistic_name = "dCov",
      statistic_function = dcov_stat,
      B = B
    )
  )
  
  counter <- counter + 1L
  
  # ----------------------------------------------------------
  # Graph-based statistic
  # ----------------------------------------------------------
  
  cat(
    "Computing graph statistic ...\n"
  )
  
  results_list[[counter]] <- cbind(
    n = n,
    analyse_one_statistic(
      dat = dat,
      statistic_name = "Graph-based",
      statistic_function = graph_stat,
      B = B
    )
  )
  
  counter <- counter + 1L
}

# ============================================================
# Combine results
# ============================================================

resampling_results <- do.call(
  rbind,
  results_list
)

rownames(
  resampling_results
) <- NULL

# ------------------------------------------------------------
# Complete results
# ------------------------------------------------------------

cat(
  "\n\n============================================\n"
)

cat(
  "BOOTSTRAP AND JACKKNIFE RESULTS\n"
)

cat(
  "============================================\n\n"
)

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

cat(
  "\n\n============================================\n"
)

cat(
  "COMPACT RESULTS\n"
)

cat(
  "============================================\n\n"
)

print(
  compact_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# Jackknife-to-bootstrap SE ratio
# ------------------------------------------------------------

compact_results$SE_Ratio <-
  compact_results$Jackknife_SE /
  compact_results$Bootstrap_SE

cat(
  "\n\n============================================\n"
)

cat(
  "STANDARD ERROR RATIO\n"
)

cat(
  "============================================\n\n"
)

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

# ============================================================
# Save results
# ============================================================

write.csv(
  resampling_results,
  file =
    "Bootstrap_Jackknife_Resampling_Results.csv",
  row.names = FALSE
)

write.csv(
  compact_results,
  file =
    "Bootstrap_Jackknife_Compact_Results.csv",
  row.names = FALSE
)

# ============================================================
# Optional bootstrap distributions
# ============================================================

bootstrap_distributions <- list()

for (n in n.grid) {
  
  cat(
    "\nGenerating bootstrap distributions for n =",
    n,
    "\n"
  )
  
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

# ============================================================
# Optional bootstrap density plots
# ============================================================

for (n in n.grid) {
  
  obj <- bootstrap_distributions[[paste0("n", n)]]
  
  par(
    mfrow = c(1, 3)
  )
  
  # tau*
  if (
    length(obj$tau_star) > 1
  ) {
    
    plot(
      density(
        obj$tau_star,
        na.rm = TRUE
      ),
      main =
        paste(
          "Bootstrap tau*",
          "\nn =",
          n
        ),
      xlab = "tau*"
    )
  }
  
  # dCov
  if (
    length(obj$dCov) > 1
  ) {
    
    plot(
      density(
        obj$dCov,
        na.rm = TRUE
      ),
      main =
        paste(
          "Bootstrap dCov",
          "\nn =",
          n
        ),
      xlab = "dCov"
    )
  }
  
  # Graph
  if (
    length(obj$graph) > 1
  ) {
    
    plot(
      density(
        obj$graph,
        na.rm = TRUE
      ),
      main =
        paste(
          "Bootstrap Graph statistic",
          "\nn =",
          n
        ),
      xlab = "Graph statistic"
    )
  }
  
  par(
    mfrow = c(1, 1)
  )
}

# ============================================================
# End of code
# ============================================================