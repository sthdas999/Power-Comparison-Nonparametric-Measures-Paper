############################################################
# COMPARATIVE POWER SIMULATION
#
# Tests:
#   1. Bergsma--Dassios tau*
#   2. Distance covariance
#   3. HHG-style graph-based statistic
#
# The HHG-style statistic is implemented directly in R.
# No HHG package is required.
############################################################

rm(list = ls())
gc()

set.seed(20260921)

############################################################
# 1. SIMULATION SETTINGS
############################################################

n.grid <- c(30, 50, 100)

B <- 2000
M <- 500
alpha <- 0.05

############################################################
# 2. TAU* PACKAGE
############################################################

if (!requireNamespace("TauStar", quietly = TRUE)) {
  install.packages(
    "TauStar",
    repos = "https://cloud.r-project.org"
  )
}

library(TauStar)

############################################################
# 3. BERGSMA--DASSIOS TAU*
############################################################

tau_star_scalar <- function(x, y) {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  ok <- is.finite(x) & is.finite(y)
  
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 5) {
    return(0)
  }
  
  z <- tryCatch(
    TauStar::tStar(
      x,
      y,
      method = "heller"
    ),
    error = function(e) {
      NA_real_
    }
  )
  
  if (!is.finite(z)) {
    return(0)
  }
  
  abs(as.numeric(z))
}

############################################################
# 4. DISTANCE COVARIANCE
############################################################

distance_covariance <- function(X, Y) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  n <- nrow(X)
  
  if (nrow(Y) != n) {
    stop("X and Y must have the same number of observations.")
  }
  
  DX <- as.matrix(dist(X))
  DY <- as.matrix(dist(Y))
  
  A <- DX -
    matrix(rowMeans(DX), n, n, byrow = TRUE) -
    matrix(colMeans(DX), n, n, byrow = FALSE) +
    mean(DX)
  
  B <- DY -
    matrix(rowMeans(DY), n, n, byrow = TRUE) -
    matrix(colMeans(DY), n, n, byrow = FALSE) +
    mean(DY)
  
  dcov2 <- mean(A * B)
  
  sqrt(max(dcov2, 0))
}

############################################################
# 5. HHG-STYLE GRAPH STATISTIC
#
# For each observation i, observations are ranked according
# to their distances from i in X-space and Y-space.
#
# For every pair of neighbourhood sizes (kx, ky), a 2 x 2
# contingency table is formed:
#
#                 Y-neighbour   Y-nonneighbour
# X-neighbour          a              b
# X-nonneighbour       c              d
#
# The local Pearson chi-square statistics are calculated,
# and the maximum is used as the global HHG-style statistic.
############################################################

hhg_statistic <- function(X, Y) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  n <- nrow(X)
  
  if (nrow(Y) != n) {
    stop("X and Y must have the same number of observations.")
  }
  
  if (n < 5) {
    return(0)
  }
  
  DX <- as.matrix(dist(X))
  DY <- as.matrix(dist(Y))
  
  max_chisq <- 0
  
  ##########################################################
  # Loop over anchor observations
  ##########################################################
  
  for (i in seq_len(n)) {
    
    ########################################################
    # Rank distances from observation i
    ########################################################
    
    rx <- rank(
      DX[i, ],
      ties.method = "first"
    )
    
    ry <- rank(
      DY[i, ],
      ties.method = "first"
    )
    
    ########################################################
    # The anchor itself must never belong to a neighbourhood
    ########################################################
    
    rx[i] <- Inf
    ry[i] <- Inf
    
    ########################################################
    # Neighbourhood sizes
    ########################################################
    
    for (kx in 1:(n - 2)) {
      
      in_x <- rx <= kx
      
      in_x[i] <- FALSE
      
      nx <- sum(in_x)
      
      if (nx == 0 || nx >= n - 1) {
        next
      }
      
      for (ky in 1:(n - 2)) {
        
        in_y <- ry <= ky
        
        in_y[i] <- FALSE
        
        ny <- sum(in_y)
        
        if (ny == 0 || ny >= n - 1) {
          next
        }
        
        ####################################################
        # 2 x 2 contingency table
        ####################################################
        
        a <- sum(in_x & in_y)
        
        b <- sum(in_x & !in_y)
        
        c <- sum(!in_x & in_y)
        
        d <- sum(!in_x & !in_y)
        
        ####################################################
        # Pearson chi-square statistic
        ####################################################
        
        denominator <-
          (a + b) *
          (c + d) *
          (a + c) *
          (b + d)
        
        if (denominator <= 0) {
          next
        }
        
        chisq <-
          (n - 1) *
          (a * d - b * c)^2 /
          denominator
        
        if (is.finite(chisq) &&
            chisq > max_chisq) {
          
          max_chisq <- chisq
        }
      }
    }
  }
  
  as.numeric(max_chisq)
}

############################################################
# 6. COORDINATEWISE MAXIMUM TAU*
#
# Used for the five-dimensional Model V.
############################################################

tau_star_max <- function(X, Y) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  p <- ncol(X)
  q <- ncol(Y)
  
  values <- numeric(p * q)
  
  z <- 1
  
  for (j in seq_len(p)) {
    
    for (k in seq_len(q)) {
      
      values[z] <-
        tau_star_scalar(
          X[, j],
          Y[, k]
        )
      
      z <- z + 1
    }
  }
  
  max(values, na.rm = TRUE)
}

############################################################
# 7. DATA-GENERATING MODELS
############################################################

generate_model <- function(model, n) {
  
  ##########################################################
  # MODEL I: LINEAR
  ##########################################################
  
  if (model == 1) {
    
    X <- rnorm(n)
    
    epsilon <- rnorm(
      n,
      mean = 0,
      sd = 0.5
    )
    
    Y <- X + epsilon
    
    return(
      list(
        X = matrix(X, ncol = 1),
        Y = matrix(Y, ncol = 1)
      )
    )
  }
  
  ##########################################################
  # MODEL II: QUADRATIC
  ##########################################################
  
  if (model == 2) {
    
    X <- rnorm(n)
    
    epsilon <- rnorm(
      n,
      mean = 0,
      sd = 0.5
    )
    
    Y <- X^2 + epsilon
    
    return(
      list(
        X = matrix(X, ncol = 1),
        Y = matrix(Y, ncol = 1)
      )
    )
  }
  
  ##########################################################
  # MODEL III: CIRCULAR
  ##########################################################
  
  if (model == 3) {
    
    theta <- runif(
      n,
      min = 0,
      max = 2 * pi
    )
    
    epsilon <- rnorm(n)
    
    R <- 1 + 0.1 * epsilon
    
    X <- R * cos(theta)
    
    Y <- R * sin(theta)
    
    return(
      list(
        X = matrix(X, ncol = 1),
        Y = matrix(Y, ncol = 1)
      )
    )
  }
  
  ##########################################################
  # MODEL IV: ORDINAL
  ##########################################################
  
  if (model == 4) {
    
    X <- rnorm(n)
    
    epsilon <- rnorm(n)
    
    Y.star <- X + epsilon
    
    Y <- ifelse(
      Y.star <= -0.5,
      1,
      ifelse(
        Y.star <= 0.5,
        2,
        3
      )
    )
    
    return(
      list(
        X = matrix(X, ncol = 1),
        Y = matrix(Y, ncol = 1)
      )
    )
  }
  
  ##########################################################
  # MODEL V: FIVE-DIMENSIONAL GAUSSIAN
  ##########################################################
  
  if (model == 5) {
    
    p <- 5
    
    rho <- 0.5
    
    X <- matrix(
      rnorm(n * p),
      nrow = n,
      ncol = p
    )
    
    epsilon <- matrix(
      rnorm(n * p),
      nrow = n,
      ncol = p
    )
    
    Y <-
      rho * X +
      sqrt(1 - rho^2) * epsilon
    
    return(
      list(
        X = X,
        Y = Y
      )
    )
  }
  
  stop("Invalid model number.")
}

############################################################
# 8. OBSERVED TEST STATISTICS
############################################################

calculate_statistics <- function(X, Y) {
  
  ##########################################################
  # TAU*
  ##########################################################
  
  if (
    ncol(X) == 1 &&
    ncol(Y) == 1
  ) {
    
    tau_obs <-
      tau_star_scalar(
        X[, 1],
        Y[, 1]
      )
    
  } else {
    
    tau_obs <-
      tau_star_max(
        X,
        Y
      )
  }
  
  ##########################################################
  # DISTANCE COVARIANCE
  ##########################################################
  
  dcov_obs <-
    distance_covariance(
      X,
      Y
    )
  
  ##########################################################
  # HHG
  ##########################################################
  
  hhg_obs <-
    hhg_statistic(
      X,
      Y
    )
  
  c(
    tau_star = tau_obs,
    dCov = dcov_obs,
    HHG = hhg_obs
  )
}

############################################################
# 9. PERMUTATION TEST
############################################################

permutation_test <- function(
    X,
    Y,
    M = 500,
    alpha = 0.05) {
  
  ##########################################################
  # Observed statistics
  ##########################################################
  
  observed <-
    calculate_statistics(
      X,
      Y
    )
  
  ##########################################################
  # Storage
  ##########################################################
  
  tau_perm <- numeric(M)
  
  dcov_perm <- numeric(M)
  
  hhg_perm <- numeric(M)
  
  ##########################################################
  # Permutations
  ##########################################################
  
  for (m in seq_len(M)) {
    
    index <- sample.int(
      nrow(Y),
      size = nrow(Y),
      replace = FALSE
    )
    
    Y_perm <-
      Y[index, , drop = FALSE]
    
    ########################################################
    # tau*
    ########################################################
    
    if (
      ncol(X) == 1 &&
      ncol(Y) == 1
    ) {
      
      tau_perm[m] <-
        tau_star_scalar(
          X[, 1],
          Y_perm[, 1]
        )
      
    } else {
      
      tau_perm[m] <-
        tau_star_max(
          X,
          Y_perm
        )
    }
    
    ########################################################
    # dCov
    ########################################################
    
    dcov_perm[m] <-
      distance_covariance(
        X,
        Y_perm
      )
    
    ########################################################
    # HHG
    ########################################################
    
    hhg_perm[m] <-
      hhg_statistic(
        X,
        Y_perm
      )
  }
  
  ##########################################################
  # Permutation p-values
  ##########################################################
  
  p_tau <-
    (
      1 +
        sum(
          tau_perm >= observed["tau_star"]
        )
    ) /
    (M + 1)
  
  p_dcov <-
    (
      1 +
        sum(
          dcov_perm >= observed["dCov"]
        )
    ) /
    (M + 1)
  
  p_hhg <-
    (
      1 +
        sum(
          hhg_perm >= observed["HHG"]
        )
    ) /
    (M + 1)
  
  ##########################################################
  # Rejection indicators
  ##########################################################
  
  reject_tau <-
    as.numeric(
      p_tau < alpha
    )
  
  reject_dcov <-
    as.numeric(
      p_dcov < alpha
    )
  
  reject_hhg <-
    as.numeric(
      p_hhg < alpha
    )
  
  c(
    p_tau_star = p_tau,
    p_dCov = p_dcov,
    p_HHG = p_hhg,
    
    reject_tau_star = reject_tau,
    reject_dCov = reject_dcov,
    reject_HHG = reject_hhg
  )
}

############################################################
# 10. ONE MONTE CARLO REPLICATION
############################################################

one_replication <- function(
    model,
    n,
    M = 500,
    alpha = 0.05) {
  
  dat <-
    generate_model(
      model = model,
      n = n
    )
  
  permutation_test(
    X = dat$X,
    Y = dat$Y,
    M = M,
    alpha = alpha
  )
}

############################################################
# 11. RUN ONE MODEL / SAMPLE SIZE
############################################################

run_condition <- function(
    model,
    n,
    B = 2000,
    M = 500,
    alpha = 0.05) {
  
  rejection_matrix <-
    matrix(
      0,
      nrow = B,
      ncol = 3
    )
  
  colnames(rejection_matrix) <-
    c(
      "tau_star",
      "dCov",
      "HHG"
    )
  
  for (b in seq_len(B)) {
    
    if (b %% 100 == 0) {
      
      cat(
        "Model =",
        model,
        "| n =",
        n,
        "| Replication =",
        b,
        "/",
        B,
        "\n"
      )
    }
    
    z <-
      one_replication(
        model = model,
        n = n,
        M = M,
        alpha = alpha
      )
    
    rejection_matrix[b, ] <-
      c(
        z["reject_tau_star"],
        z["reject_dCov"],
        z["reject_HHG"]
      )
  }
  
  ##########################################################
  # Empirical power
  ##########################################################
  
  power <-
    colMeans(
      rejection_matrix
    )
  
  ##########################################################
  # Monte Carlo standard error
  ##########################################################
  
  MCSE <-
    sqrt(
      power * (1 - power) / B
    )
  
  data.frame(
    Model = model,
    n = n,
    
    tau_star = unname(
      power["tau_star"]
    ),
    
    dCov = unname(
      power["dCov"]
    ),
    
    HHG = unname(
      power["HHG"]
    ),
    
    MCSE_tau_star = unname(
      MCSE["tau_star"]
    ),
    
    MCSE_dCov = unname(
      MCSE["dCov"]
    ),
    
    MCSE_HHG = unname(
      MCSE["HHG"]
    ),
    
    row.names = NULL
  )
}

############################################################
# 12. COMPLETE SIMULATION
############################################################

results_list <- list()

counter <- 1

for (model in 1:5) {
  
  for (n in n.grid) {
    
    cat(
      "\n============================================\n"
    )
    
    cat(
      "Model:",
      model,
      "| Sample size:",
      n,
      "\n"
    )
    
    cat(
      "============================================\n"
    )
    
    results_list[[counter]] <-
      run_condition(
        model = model,
        n = n,
        B = B,
        M = M,
        alpha = alpha
      )
    
    counter <- counter + 1
  }
}

############################################################
# 13. COMBINE RESULTS
############################################################

results <-
  do.call(
    rbind,
    results_list
  )

rownames(results) <- NULL

############################################################
# 14. MODEL LABELS
############################################################

results$Model_Name <-
  c(
    "Linear",
    "Linear",
    "Linear",
    "Quadratic",
    "Quadratic",
    "Quadratic",
    "Circular",
    "Circular",
    "Circular",
    "Ordinal",
    "Ordinal",
    "Ordinal",
    "5D Gaussian",
    "5D Gaussian",
    "5D Gaussian"
  )

############################################################
# 15. REORDER RESULTS
############################################################

results <-
  results[
    ,
    c(
      "Model",
      "Model_Name",
      "n",
      "tau_star",
      "dCov",
      "HHG",
      "MCSE_tau_star",
      "MCSE_dCov",
      "MCSE_HHG"
    )
  ]

############################################################
# 16. DISPLAY RESULTS
############################################################

print(
  results,
  row.names = FALSE
)

############################################################
# 17. SAVE COMPLETE RESULTS
############################################################

write.csv(
  results,
  "comparative_power_results_tau_dCov_HHG.csv",
  row.names = FALSE
)

############################################################
# 18. COMPACT POWER TABLE
############################################################

power_table <-
  results[
    ,
    c(
      "Model_Name",
      "n",
      "tau_star",
      "dCov",
      "HHG"
    )
  ]

print(
  power_table,
  row.names = FALSE
)

############################################################
# END
############################################################
