# ============================================================
# Bootstrap Power Analysis
# ============================================================

rm(list = ls())

set.seed(2026)

# ------------------------------------------------------------
# Required package
# ------------------------------------------------------------
library(energy)

# ------------------------------------------------------------
# Simulation settings
# ------------------------------------------------------------
alpha   <- 0.05
B_boot  <- 2000
N_sim   <- 1000
n_grid  <- c(50, 100, 200)

# ============================================================
# 1. Bergsma--Dassios tau* statistic
# ============================================================

tau_star_stat <- function(x, y) {
  
  n <- length(x)
  
  # Convert observations to ranks
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  
  z <- 0
  
  # Fourth-order summation
  for (i in 1:(n - 3)) {
    
    for (j in (i + 1):(n - 2)) {
      
      for (k in (j + 1):(n - 1)) {
        
        for (l in (k + 1):n) {
          
          # Four observations
          a <- c(rx[i], rx[j], rx[k], rx[l])
          b <- c(ry[i], ry[j], ry[k], ry[l])
          
          # Ordering according to X and Y
          ox <- order(a)
          oy <- order(b)
          
          x1 <- a[ox[1]]
          x2 <- a[ox[2]]
          x3 <- a[ox[3]]
          x4 <- a[ox[4]]
          
          y1 <- b[oy[1]]
          y2 <- b[oy[2]]
          y3 <- b[oy[3]]
          y4 <- b[oy[4]]
          
          # Pairwise concordance/discordance terms
          s1 <- sign((x1 - x2) * (y1 - y2))
          s2 <- sign((x3 - x4) * (y3 - y4))
          
          s3 <- sign((x1 - x3) * (y1 - y3))
          s4 <- sign((x2 - x4) * (y2 - y4))
          
          s5 <- sign((x1 - x4) * (y1 - y4))
          s6 <- sign((x2 - x3) * (y2 - y3))
          
          z <- z +
            (s1 * s2 +
               s3 * s4 +
               s5 * s6) / 3
        }
      }
    }
  }
  
  # Normalization
  z / choose(n, 4)
}


# ============================================================
# 2. Distance covariance statistic
# ============================================================

dcov_stat <- function(x, y) {
  
  n <- length(x)
  
  # Pairwise Euclidean distances
  dx <- as.matrix(dist(x))
  dy <- as.matrix(dist(y))
  
  # Double centering for X
  Ax <- dx -
    matrix(rowMeans(dx), n, n) -
    matrix(colMeans(dx), n, n, byrow = TRUE) +
    mean(dx)
  
  # Double centering for Y
  Ay <- dy -
    matrix(rowMeans(dy), n, n) -
    matrix(colMeans(dy), n, n, byrow = TRUE) +
    mean(dy)
  
  # Sample distance covariance
  sqrt(mean(Ax * Ay))
}


# ============================================================
# 3. Graph-based nearest-neighbour statistic
# ============================================================

graph_stat <- function(x, y) {
  
  n <- length(x)
  
  # Distance matrix based on X
  D <- as.matrix(dist(x))
  
  # Exclude self-neighbours
  diag(D) <- Inf
  
  # Nearest neighbour of each observation
  nn <- apply(D, 1, which.min)
  
  statistic <- numeric(n)
  
  # Local concordance indicator
  for (i in 1:n) {
    
    j <- nn[i]
    
    statistic[i] <-
      as.numeric(
        sign(x[i] - x[j]) *
          sign(y[i] - y[j]) > 0
      )
  }
  
  # Average local concordance
  mean(statistic)
}


# ============================================================
# 4. Bootstrap / permutation null distribution
# ============================================================

bootstrap_critical <- function(x,
                               y,
                               statistic,
                               B = 2000,
                               alpha = 0.05) {
  
  n <- length(x)
  
  null_stat <- numeric(B)
  
  for (b in 1:B) {
    
    # Permute Y relative to X
    y_perm <- sample(
      y,
      size = n,
      replace = FALSE
    )
    
    # Calculate statistic under the null
    null_stat[b] <-
      statistic(x, y_perm)
  }
  
  # Upper (1-alpha) critical value
  quantile(
    null_stat,
    probs = 1 - alpha,
    names = FALSE,
    na.rm = TRUE
  )
}


# ============================================================
# 5. Bootstrap power for a given sample size and statistic
# ============================================================

bootstrap_power <- function(n,
                            statistic,
                            B = 2000,
                            N = 1000,
                            alpha = 0.05) {
  
  # Store rejection indicators
  reject <- numeric(N)
  
  # ----------------------------------------------------------
  # Monte Carlo repetitions
  # ----------------------------------------------------------
  
  for (r in 1:N) {
    
    # --------------------------------------------------------
    # Generate data under the nonlinear alternative
    #
    # Y = X^2 + epsilon
    #
    # X ~ N(0,1)
    # epsilon ~ N(0,0.5^2)
    # --------------------------------------------------------
    
    x <- rnorm(
      n,
      mean = 0,
      sd = 1
    )
    
    epsilon <- rnorm(
      n,
      mean = 0,
      sd = 0.5
    )
    
    y <- x^2 + epsilon
    
    
    # --------------------------------------------------------
    # Observed test statistic
    # --------------------------------------------------------
    
    T_obs <- statistic(x, y)
    
    
    # --------------------------------------------------------
    # Generate null distribution by permutation
    # --------------------------------------------------------
    
    T_null <- numeric(B)
    
    for (b in 1:B) {
      
      y_perm <- sample(
        y,
        size = n,
        replace = FALSE
      )
      
      T_null[b] <-
        statistic(x, y_perm)
    }
    
    
    # --------------------------------------------------------
    # Bootstrap/permutation critical value
    # --------------------------------------------------------
    
    critical_value <- quantile(
      T_null,
      probs = 1 - alpha,
      names = FALSE,
      na.rm = TRUE
    )
    
    
    # --------------------------------------------------------
    # Rejection indicator
    # --------------------------------------------------------
    
    reject[r] <-
      as.numeric(
        T_obs > critical_value
      )
  }
  
  # ----------------------------------------------------------
  # Estimated power
  # ----------------------------------------------------------
  
  mean(reject)
}


# ============================================================
# 6. Run the complete bootstrap power simulation
# ============================================================

results <- data.frame(
  n = integer(),
  Method = character(),
  Power = numeric(),
  stringsAsFactors = FALSE
)


for (n in n_grid) {
  
  cat("\n============================================\n")
  cat("Sample size:", n, "\n")
  cat("============================================\n")
  
  
  # ----------------------------------------------------------
  # Bergsma--Dassios tau*
  # ----------------------------------------------------------
  
  cat("Running tau* ...\n")
  
  p_tau <- bootstrap_power(
    n = n,
    statistic = tau_star_stat,
    B = B_boot,
    N = N_sim,
    alpha = alpha
  )
  
  results <- rbind(
    results,
    data.frame(
      n = n,
      Method = "Tau-star",
      Power = p_tau
    )
  )
  
  
  # ----------------------------------------------------------
  # Distance covariance
  # ----------------------------------------------------------
  
  cat("Running distance covariance ...\n")
  
  p_dcov <- bootstrap_power(
    n = n,
    statistic = dcov_stat,
    B = B_boot,
    N = N_sim,
    alpha = alpha
  )
  
  results <- rbind(
    results,
    data.frame(
      n = n,
      Method = "dCov",
      Power = p_dcov
    )
  )
  
  
  # ----------------------------------------------------------
  # Graph-based statistic
  # ----------------------------------------------------------
  
  cat("Running graph-based statistic ...\n")
  
  p_graph <- bootstrap_power(
    n = n,
    statistic = graph_stat,
    B = B_boot,
    N = N_sim,
    alpha = alpha
  )
  
  results <- rbind(
    results,
    data.frame(
      n = n,
      Method = "Graph-based",
      Power = p_graph
    )
  )
}


# ============================================================
# 7. Round and display power estimates
# ============================================================

results$Power <-
  round(
    results$Power,
    4
  )

print(results)


# ============================================================
# 8. Convert results to wide format
# ============================================================

power_table <- reshape(
  results,
  idvar = "n",
  timevar = "Method",
  direction = "wide"
)

# Remove the "Power." prefix
names(power_table) <-
  sub(
    "^Power\\.",
    "",
    names(power_table)
  )

print(power_table)


# ============================================================
# 9. Save results
# ============================================================

write.csv(
  results,
  "bootstrap_power_results.csv",
  row.names = FALSE
)

write.csv(
  power_table,
  "bootstrap_power_table.csv",
  row.names = FALSE
)


# ============================================================
# 10. Optional compact display
# ============================================================

cat("\n============================================\n")
cat("Bootstrap Power Results\n")
cat("============================================\n")

print(
  power_table,
  row.names = FALSE
)