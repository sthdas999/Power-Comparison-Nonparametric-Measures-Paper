############################################################
# COMPARATIVE POWER STUDY
#
# Tests:
#   1. Bergsma--Dassios tau*
#   2. Distance covariance
#   3. Heller--Heller--Gorfine (HHG)
#
# Models:
#   I.   Linear dependence
#   II.  Quadratic dependence
#   III. Circular dependence
#   IV.  Ordinal dependence
#   V.   Five-dimensional Gaussian dependence
#
# Sample sizes:
#   n = 30, 50, 100
#
# Monte Carlo replications:
#   B = 5000
#
# Permutations:
#   M = 1000
#
# Significance level:
#   alpha = 0.05
#
# IMPORTANT:
# This code is provided for execution by the user.
# No simulation results are generated in this manuscript.
############################################################


############################################################
# 1. INSTALL AND LOAD REQUIRED PACKAGES
############################################################

# TauStar
if (!requireNamespace("TauStar", quietly = TRUE)) {
  install.packages("TauStar")
}

# HHG
#
# HHG has been archived from the current CRAN repository.
# Therefore, install the archived CRAN version if necessary.
############################################################

if (!requireNamespace("HHG", quietly = TRUE)) {
  
  hhg_url <-
    "https://cran.r-project.org/src/contrib/Archive/HHG/HHG_2.3.7.tar.gz"
  
  install.packages(
    hhg_url,
    repos = NULL,
    type = "source"
  )
}

library(TauStar)
library(HHG)


############################################################
# 2. SIMULATION SETTINGS
############################################################

set.seed(20260921)

n.grid <- c(30, 50, 100)

B <- 5000

M <- 1000

alpha <- 0.05


############################################################
# 3. MODEL I
#    Linear dependence
#
#    Y = X + epsilon
#
#    X ~ N(0,1)
#    epsilon ~ N(0,0.5^2)
############################################################

generate_model_1 <- function(n) {
  
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
  
  Y <- X + epsilon
  
  return(
    list(
      X = X,
      Y = Y
    )
  )
}


############################################################
# 4. MODEL II
#    Quadratic dependence
#
#    Y = X^2 + epsilon
#
#    X ~ N(0,1)
#    epsilon ~ N(0,0.5^2)
############################################################

generate_model_2 <- function(n) {
  
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
  
  return(
    list(
      X = X,
      Y = Y
    )
  )
}


############################################################
# 5. MODEL III
#    Circular dependence
#
#    X = R cos(theta)
#    Y = R sin(theta)
#
#    R = 1 + 0.1 epsilon
#    theta ~ Uniform(0,2*pi)
############################################################

generate_model_3 <- function(n) {
  
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
      X = X,
      Y = Y
    )
  )
}


############################################################
# 6. MODEL IV
#    Ordinal dependence
#
#    Y* = X + epsilon
#
#    Y = 1 if Y* <= -0.5
#        2 if -0.5 < Y* <= 0.5
#        3 if Y* > 0.5
############################################################

generate_model_4 <- function(n) {
  
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
      X = X,
      Y = Y
    )
  )
}


############################################################
# 7. MODEL V
#    Five-dimensional Gaussian dependence
#
#    X ~ N_5(0,I_5)
#    epsilon ~ N_5(0,I_5)
#
#    Y = rho X + sqrt(1-rho^2) epsilon
#
#    rho = 0.5
############################################################

generate_model_5 <- function(
    n,
    rho = 0.5) {
  
  X <- matrix(
    rnorm(5 * n),
    nrow = n,
    ncol = 5
  )
  
  epsilon <- matrix(
    rnorm(5 * n),
    nrow = n,
    ncol = 5
  )
  
  Y <- rho * X +
    sqrt(1 - rho^2) * epsilon
  
  colnames(X) <- paste0(
    "X",
    1:5
  )
  
  colnames(Y) <- paste0(
    "Y",
    1:5
  )
  
  return(
    list(
      X = X,
      Y = Y
    )
  )
}


############################################################
# 8. GENERAL DATA-GENERATING FUNCTION
############################################################

generate_data <- function(
    model,
    n) {
  
  if (model == 1) {
    
    return(
      generate_model_1(n)
    )
    
  } else if (model == 2) {
    
    return(
      generate_model_2(n)
    )
    
  } else if (model == 3) {
    
    return(
      generate_model_3(n)
    )
    
  } else if (model == 4) {
    
    return(
      generate_model_4(n)
    )
    
  } else if (model == 5) {
    
    return(
      generate_model_5(n)
    )
    
  } else {
    
    stop(
      "Model must be one of 1, 2, 3, 4, or 5."
    )
  }
}


############################################################
# 9. MATRIX CONVERSION FUNCTION
############################################################

as_data_matrix <- function(x) {
  
  if (is.null(dim(x))) {
    
    return(
      matrix(
        x,
        ncol = 1
      )
    )
    
  } else {
    
    return(
      as.matrix(x)
    )
  }
}


############################################################
# 10. BERGSMA--DASSIOS TAU* STATISTIC
#
# For Models I--IV:
#   standard scalar tau*
#
# For Model V:
#   maximum absolute coordinatewise tau*
############################################################

tau_star_scalar <- function(
    x,
    y) {
  
  as.numeric(
    TauStar::tStar(
      x,
      y,
      method = "heller"
    )
  )
}


############################################################
# 11. COORDINATEWISE MAXIMUM TAU*
#
# Used only for the five-dimensional model.
#
# T_tau,max =
# max_{j,k} |tau*(X_j,Y_k)|
############################################################

tau_star_max <- function(
    X,
    Y) {
  
  X <- as_data_matrix(X)
  
  Y <- as_data_matrix(Y)
  
  p <- ncol(X)
  
  q <- ncol(Y)
  
  tau.matrix <- matrix(
    NA_real_,
    nrow = p,
    ncol = q
  )
  
  for (j in seq_len(p)) {
    
    for (k in seq_len(q)) {
      
      tau.matrix[j, k] <-
        tau_star_scalar(
          X[, j],
          Y[, k]
        )
    }
  }
  
  return(
    max(
      abs(tau.matrix)
    )
  )
}


############################################################
# 12. TAU* OBSERVED STATISTIC
############################################################

tau_statistic <- function(
    X,
    Y) {
  
  X <- as_data_matrix(X)
  
  Y <- as_data_matrix(Y)
  
  if (
    ncol(X) == 1 &&
    ncol(Y) == 1
  ) {
    
    return(
      abs(
        tau_star_scalar(
          X[, 1],
          Y[, 1]
        )
      )
    )
    
  } else {
    
    return(
      tau_star_max(
        X,
        Y
      )
    )
  }
}


############################################################
# 13. DISTANCE COVARIANCE STATISTIC
#
# The distance covariance is computed directly from the
# doubly centered Euclidean distance matrices.
#
# No additional package is required.
############################################################

dcov_statistic <- function(
    X,
    Y) {
  
  X <- as_data_matrix(X)
  
  Y <- as_data_matrix(Y)
  
  n <- nrow(X)
  
  
  ##########################################################
  # Euclidean distance matrices
  ##########################################################
  
  A <- as.matrix(
    dist(X)
  )
  
  B <- as.matrix(
    dist(Y)
  )
  
  
  ##########################################################
  # Double centering
  ##########################################################
  
  A <- A -
    matrix(
      rowMeans(A),
      nrow = n,
      ncol = n
    ) -
    matrix(
      colMeans(A),
      nrow = n,
      ncol = n,
      byrow = TRUE
    ) +
    mean(A)
  
  
  B <- B -
    matrix(
      rowMeans(B),
      nrow = n,
      ncol = n
    ) -
    matrix(
      colMeans(B),
      nrow = n,
      ncol = n,
      byrow = TRUE
    ) +
    mean(B)
  
  
  ##########################################################
  # Sample distance covariance
  ##########################################################
  
  dcov2 <- mean(
    A * B
  )
  
  
  return(
    sqrt(
      max(
        dcov2,
        0
      )
    )
  )
}


############################################################
# 14. DISTANCE MATRICES FOR HHG
############################################################

distance_matrix <- function(X) {
  
  X <- as_data_matrix(X)
  
  as.matrix(
    dist(X)
  )
}


############################################################
# 15. HHG TEST
#
# HHG accepts distance matrices as its input.
#
# The statistic used here is the sum of Pearson
# chi-square statistics, which is one of the standard
# HHG statistics.
############################################################

hhg_test <- function(
    X,
    Y,
    M = 1000) {
  
  X <- as_data_matrix(X)
  
  Y <- as_data_matrix(Y)
  
  Dx <- distance_matrix(X)
  
  Dy <- distance_matrix(Y)
  
  result <- HHG::hhg.test(
    Dx,
    Dy,
    ties = TRUE,
    w.sum = 0,
    w.max = 2,
    nr.perm = M,
    is.sequential = FALSE,
    tables.wanted = FALSE,
    perm.stats.wanted = FALSE
  )
  
  return(
    as.numeric(
      result$perm.pval.hhg.sc
    )
  )
}


############################################################
# 16. GENERIC PERMUTATION TEST
#
# Used for tau* and distance covariance.
############################################################

permutation_test <- function(
    X,
    Y,
    statistic_function,
    M = 1000) {
  
  X <- as_data_matrix(X)
  
  Y <- as_data_matrix(Y)
  
  n <- nrow(X)
  
  observed <- statistic_function(
    X,
    Y
  )
  
  permuted_statistics <- numeric(M)
  
  for (m in seq_len(M)) {
    
    permutation <- sample.int(
      n,
      size = n,
      replace = FALSE
    )
    
    Y.permuted <-
      Y[permutation, , drop = FALSE]
    
    permuted_statistics[m] <-
      statistic_function(
        X,
        Y.permuted
      )
  }
  
  p.value <- (
    1 +
      sum(
        permuted_statistics >=
          observed
      )
  ) /
    (M + 1)
  
  return(
    p.value
  )
}


############################################################
# 17. TAU* PERMUTATION TEST
############################################################

tau_star_test <- function(
    X,
    Y,
    M = 1000) {
  
  permutation_test(
    X = X,
    Y = Y,
    statistic_function = tau_statistic,
    M = M
  )
}


############################################################
# 18. DISTANCE COVARIANCE PERMUTATION TEST
############################################################

dcov_test <- function(
    X,
    Y,
    M = 1000) {
  
  permutation_test(
    X = X,
    Y = Y,
    statistic_function = dcov_statistic,
    M = M
  )
}


############################################################
# 19. ONE MONTE CARLO REPLICATION
############################################################

one_replication <- function(
    model,
    n,
    M = 1000,
    alpha = 0.05) {
  
  ##########################################################
  # Generate data
  ##########################################################
  
  data <- generate_data(
    model = model,
    n = n
  )
  
  X <- data$X
  
  Y <- data$Y
  
  
  ##########################################################
  # Bergsma--Dassios tau*
  ##########################################################
  
  p_tau <- tau_star_test(
    X = X,
    Y = Y,
    M = M
  )
  
  
  ##########################################################
  # Distance covariance
  ##########################################################
  
  p_dcov <- dcov_test(
    X = X,
    Y = Y,
    M = M
  )
  
  
  ##########################################################
  # HHG
  ##########################################################
  
  p_hhg <- hhg_test(
    X = X,
    Y = Y,
    M = M
  )
  
  
  ##########################################################
  # Rejection indicators
  ##########################################################
  
  rejection_tau <-
    as.numeric(
      p_tau <= alpha
    )
  
  rejection_dcov <-
    as.numeric(
      p_dcov <= alpha
    )
  
  rejection_hhg <-
    as.numeric(
      p_hhg <= alpha
    )
  
  
  ##########################################################
  # Return results
  ##########################################################
  
  return(
    c(
      tau_star = rejection_tau,
      dCov = rejection_dcov,
      HHG = rejection_hhg
    )
  )
}


############################################################
# 20. RUN ONE MODEL AND ONE SAMPLE SIZE
############################################################

run_condition <- function(
    model,
    n,
    B = 5000,
    M = 1000,
    alpha = 0.05) {
  
  rejection_matrix <- matrix(
    0,
    nrow = B,
    ncol = 3
  )
  
  colnames(
    rejection_matrix
  ) <- c(
    "tau_star",
    "dCov",
    "HHG"
  )
  
  
  for (b in seq_len(B)) {
    
    rejection_matrix[b, ] <-
      one_replication(
        model = model,
        n = n,
        M = M,
        alpha = alpha
      )
    
    
    if (
      b %% 100 == 0
    ) {
      
      cat(
        "Model:",
        model,
        "| n:",
        n,
        "| Replication:",
        b,
        "of",
        B,
        "\n"
      )
    }
  }
  
  
  ##########################################################
  # Empirical power
  ##########################################################
  
  power <- colMeans(
    rejection_matrix
  )
  
  
  return(
    data.frame(
      Model = model,
      Sample_Size = n,
      tau_star = power["tau_star"],
      dCov = power["dCov"],
      HHG = power["HHG"],
      row.names = NULL
    )
  )
}


############################################################
# 21. COMPLETE SIMULATION
############################################################

run_complete_simulation <- function(
    n.grid = c(30, 50, 100),
    B = 5000,
    M = 1000,
    alpha = 0.05) {
  
  all_results <- list()
  
  counter <- 1
  
  
  for (model in 1:5) {
    
    for (n in n.grid) {
      
      cat(
        "\n====================================\n"
      )
      
      cat(
        "Starting Model:",
        model,
        "| Sample size:",
        n,
        "\n"
      )
      
      cat(
        "====================================\n"
      )
      
      
      all_results[[counter]] <-
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
  
  
  final_results <- do.call(
    rbind,
    all_results
  )
  
  
  rownames(
    final_results
  ) <- NULL
  
  
  return(
    final_results
  )
}


############################################################
# 22. RUN THE SIMULATION
#
# IMPORTANT:
# This command is intentionally commented out.
#
# Uncomment it only when you want to execute the
# complete simulation.
############################################################

simulation_results <-
  run_complete_simulation(
    n.grid = c(30, 50, 100),
    B = 5000,
    M = 1000,
    alpha = 0.05
  )


############################################################
# 23. SAVE THE RESULTS
#
# Execute after the simulation has been completed.
############################################################

write.csv(
  simulation_results,
  file = "comparative_power_results.csv",
  row.names = FALSE
)


############################################################
# 24. OPTIONAL: DISPLAY THE RESULTS
############################################################

print(
  simulation_results
)
