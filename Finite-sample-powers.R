## ============================================================
## FINITE-SAMPLE POWER STUDY
## Model: Y = X + epsilon
## X ~ N(0,1), epsilon ~ N(0,0.5^2), X independent of epsilon
## ============================================================

rm(list = ls())

set.seed(12345)

## ------------------------------------------------------------
## Simulation settings
## ------------------------------------------------------------

n.grid <- c(
  50, 100, 200, 500, 800,
  1000, 2000, 5000, 10000
)

B <- 1000
M <- 500
alpha <- 0.05


## ============================================================
## 1. DATA GENERATING MECHANISM
## ============================================================

generate.data <- function(n) {
  
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
  
  list(
    X = X,
    Y = Y
  )
}


## ============================================================
## 2. BERGSMA--DASSIOS TAU-STAR
## ============================================================

tau.star.stat <- function(X, Y) {
  
  n <- length(X)
  
  if (length(Y) != n) {
    stop("X and Y must have the same length.")
  }
  
  if (n < 4) {
    stop("At least four observations are required.")
  }
  
  ## ----------------------------------------------------------
  ## For continuous data, observations can be ordered by X.
  ## The tau-star statistic is based on quadruples.
  ## ----------------------------------------------------------
  
  ## Obtain all quadruples
  Q <- combn(n, 4)
  
  C <- 0
  D <- 0
  
  for (j in seq_len(ncol(Q))) {
    
    id <- Q[, j]
    
    x4 <- X[id]
    y4 <- Y[id]
    
    ## Ignore quadruples containing ties
    if (
      length(unique(x4)) < 4 ||
      length(unique(y4)) < 4
    ) {
      next
    }
    
    ## Order according to X
    ord <- order(x4)
    
    y <- y4[ord]
    
    ## Three possible pairings
    s1 <- sign(
      (y[1] - y[2]) *
        (y[3] - y[4])
    )
    
    s2 <- sign(
      (y[1] - y[3]) *
        (y[2] - y[4])
    )
    
    s3 <- sign(
      (y[1] - y[4]) *
        (y[2] - y[3])
    )
    
    ## Quadruple score
    score <- s1 + s2 + s3
    
    if (score > 0) {
      C <- C + 1
    }
    
    if (score < 0) {
      D <- D + 1
    }
  }
  
  ## Normalized statistic
  tau.star <- 4 * (C - D) / choose(n, 4)
  
  return(tau.star)
}


## ============================================================
## 3. DISTANCE COVARIANCE STATISTIC
## ============================================================

dcov.stat <- function(X, Y) {
  
  n <- length(X)
  
  if (length(Y) != n) {
    stop("X and Y must have the same length.")
  }
  
  ## Pairwise Euclidean distance matrices
  A <- abs(
    outer(X, X, "-")
  )
  
  B <- abs(
    outer(Y, Y, "-")
  )
  
  ## Double centering
  A <- sweep(
    A,
    1,
    rowMeans(A),
    "-"
  )
  
  A <- sweep(
    A,
    2,
    colMeans(A),
    "-"
  )
  
  A <- A + mean(A)
  
  B <- sweep(
    B,
    1,
    rowMeans(B),
    "-"
  )
  
  B <- sweep(
    B,
    2,
    colMeans(B),
    "-"
  )
  
  B <- B + mean(B)
  
  ## Distance covariance squared
  dcov2 <- sum(A * B) / n^2
  
  return(dcov2)
}


## ============================================================
## 4. GRAPH-BASED INDEPENDENCE STATISTIC
## ============================================================

graph.stat <- function(X, Y) {
  
  n <- length(X)
  
  if (length(Y) != n) {
    stop("X and Y must have the same length.")
  }
  
  ## Rank transformation
  RX <- rank(
    X,
    ties.method = "average"
  )
  
  RY <- rank(
    Y,
    ties.method = "average"
  )
  
  ## Pairwise distances in the bivariate rank space
  D <- as.matrix(
    dist(
      cbind(RX, RY)
    )
  )
  
  ## Construct the minimum spanning tree
  ## using Prim's algorithm
  
  selected <- rep(
    FALSE,
    n
  )
  
  selected[1] <- TRUE
  
  edges <- matrix(
    0,
    nrow = n - 1,
    ncol = 2
  )
  
  for (k in 1:(n - 1)) {
    
    best.distance <- Inf
    best.i <- NA
    best.j <- NA
    
    for (i in which(selected)) {
      
      candidates <- which(!selected)
      
      if (length(candidates) == 0) {
        break
      }
      
      j <- candidates[
        which.min(D[i, candidates])
      ]
      
      if (D[i, j] < best.distance) {
        
        best.distance <- D[i, j]
        best.i <- i
        best.j <- j
      }
    }
    
    edges[k, ] <- c(
      best.i,
      best.j
    )
    
    selected[best.j] <- TRUE
  }
  
  ## ----------------------------------------------------------
  ## Measure agreement of X and Y ordering along MST edges
  ## ----------------------------------------------------------
  
  edge.score <- numeric(
    nrow(edges)
  )
  
  for (k in seq_len(nrow(edges))) {
    
    i <- edges[k, 1]
    j <- edges[k, 2]
    
    edge.score[k] <- sign(
      (RX[i] - RX[j]) *
        (RY[i] - RY[j])
    )
  }
  
  ## Graph-based statistic
  statistic <- mean(
    edge.score > 0
  )
  
  return(statistic)
}


## ============================================================
## 5. PERMUTATION CRITICAL VALUE
## ============================================================

permutation.test <- function(
    X,
    Y,
    statistic,
    M = 500,
    alpha = 0.05
) {
  
  ## Observed statistic
  T.obs <- statistic(
    X,
    Y
  )
  
  ## Permutation statistics
  T.perm <- numeric(M)
  
  for (m in seq_len(M)) {
    
    Y.perm <- sample(
      Y,
      size = length(Y),
      replace = FALSE
    )
    
    T.perm[m] <- statistic(
      X,
      Y.perm
    )
  }
  
  ## Empirical critical value
  critical.value <- quantile(
    T.perm,
    probs = 1 - alpha,
    names = FALSE,
    type = 8
  )
  
  ## Rejection indicator
  reject <- as.integer(
    T.obs > critical.value
  )
  
  return(
    list(
      statistic = T.obs,
      critical.value = critical.value,
      reject = reject
    )
  )
}


## ============================================================
## 6. ONE MONTE CARLO REPLICATION
## ============================================================

one.replication <- function(
    n,
    M = 500,
    alpha = 0.05
) {
  
  ## Generate data
  dat <- generate.data(n)
  
  X <- dat$X
  Y <- dat$Y
  
  ## ----------------------------------------------------------
  ## Bergsma--Dassios tau-star
  ## ----------------------------------------------------------
  
  tau.result <- permutation.test(
    X = X,
    Y = Y,
    statistic = tau.star.stat,
    M = M,
    alpha = alpha
  )
  
  ## ----------------------------------------------------------
  ## Distance covariance
  ## ----------------------------------------------------------
  
  dcov.result <- permutation.test(
    X = X,
    Y = Y,
    statistic = dcov.stat,
    M = M,
    alpha = alpha
  )
  
  ## ----------------------------------------------------------
  ## Graph-based statistic
  ## ----------------------------------------------------------
  
  graph.result <- permutation.test(
    X = X,
    Y = Y,
    statistic = graph.stat,
    M = M,
    alpha = alpha
  )
  
  ## Return only rejection indicators
  c(
    tau.star = tau.result$reject,
    dCov = dcov.result$reject,
    graph = graph.result$reject
  )
}


## ============================================================
## 7. POWER FOR A GIVEN SAMPLE SIZE
## ============================================================

power.one.n <- function(
    n,
    B = 1000,
    M = 500,
    alpha = 0.05
) {
  
  rejection <- matrix(
    0,
    nrow = B,
    ncol = 3
  )
  
  colnames(rejection) <- c(
    "tau.star",
    "dCov",
    "graph"
  )
  
  for (b in seq_len(B)) {
    
    rejection[b, ] <- one.replication(
      n = n,
      M = M,
      alpha = alpha
    )
    
    if (b %% 100 == 0) {
      
      cat(
        "n =", n,
        "| replication =", b,
        "of", B,
        "\n"
      )
    }
  }
  
  ## Empirical rejection probabilities
  colMeans(
    rejection
  )
}


## ============================================================
## 8. COMPLETE SIMULATION
## ============================================================

power.results <- vector(
  "list",
  length(n.grid)
)

for (j in seq_along(n.grid)) {
  
  n <- n.grid[j]
  
  cat("\n")
  cat("============================================\n")
  cat("Sample size:", n, "\n")
  cat("============================================\n")
  
  power.results[[j]] <- power.one.n(
    n = n,
    B = B,
    M = M,
    alpha = alpha
  )
}


## ============================================================
## 9. COMBINE THE RESULTS
## ============================================================

power.table <- data.frame(
  n = n.grid,
  do.call(
    rbind,
    power.results
  ),
  row.names = NULL
)


## ============================================================
## 10. DISPLAY THE POWER TABLE
## ============================================================

print(
  power.table
)


## ============================================================
## 11. SAVE THE RESULTS
## ============================================================

write.csv(
  power.table,
  file = "finite_sample_power_results.csv",
  row.names = FALSE
)