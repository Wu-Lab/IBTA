# Inferring Co-hub Differential Networks 
source("./utils/theta2partial.R")
library(EMDomics)
library(Rfast)
update_theta_k <- function(A, alpha, ifparallel = FALSE) {
  K <- dim(A)[1]
  p <- dim(A)[2]
  theta_k <- array(rep(0, K * p * p), c(K, p, p))
  if (ifparallel) {
    tmp <- foreach(k = 1:K) %dopar% {
      results <- eigen(A[k, , ], symmetric = TRUE)
      U <- results$vectors
      D <- results$values
      D <- (-D + sqrt(D^2 + 4 * alpha[k])) / 2 / alpha[k]
      theta <- U %*% diag(D) %*% t(U)
      return(theta)
    }
    for (k in 1:K) {
      theta_k[k, , ] <- tmp[[k]]
    }
  } else {
    for (k in 1:K) {
      results <- eigen(A[k, , ], symmetric = TRUE)
      U <- results$vectors
      D <- results$values
      D <- (-D + sqrt(D^2 + 4 * alpha[k])) / 2 / alpha[k]
      theta <- U %*% diag(D) %*% t(U)
      theta_k[k, , ] <- theta
    }
  }
  return(theta_k)
}

soft_threshold <- function(X, lambda) {
  Y <- sign(X) * (pmax(abs(X) - lambda, 0))
  return(Y)
}

update_V_new <- function(C, rho, lambda2) {
  K <- dim(C)[1]
  p <- dim(C)[2]
  V <- array(rep(0, K * p * p), c(K, p, p))
  for (j in 1:p) {
    C_norm <- norm(as.vector(C[, , j]), "2")
    if (C_norm >= lambda2[j]) {
      V[, , j] <- (C_norm - lambda2[j]) / (2 * rho * C_norm) * C[, , j]
    }
  }
  return(V)
}

calculate_w <- function(X, K, p, ifparallel = FALSE) {
  w <- array(dim = c(K, p))
  if (ifparallel) {
    tmp <- foreach(
      k = 1:K,
      .combine = "cbind",
      .packages = c("EMDomics")
    ) %dopar% {
      data <- cbind(t(X[[k]]), t(X[[k + K]]))
      n1 <- dim(X[[k]])[1]
      n2 <- dim(X[[k + K]])[1]
      rownames(data) <- paste("gene", 1:p, sep = "")
      colnames(data) <- paste("sample", 1:(n1 + n2), sep = "")
      outcomes <- c(rep("A", n1), rep("B", n2))
      names(outcomes) <- colnames(data)
      results <-
        EMDomics::calculate_emd(
          data,
          outcomes,
          nperm = 10,
          parallel = FALSE,
          verbose = 0
        )
      return(results$emd[, 2])
    }
    w <- t(tmp)
  } else {
    for (k in 1:K) {
      data <- cbind(t(X[[k]]), t(X[[k + K]]))
      n1 <- dim(X[[k]])[1]
      n2 <- dim(X[[k + K]])[1]
      rownames(data) <- paste("gene", 1:p, sep = "")
      colnames(data) <- paste("sample", 1:(n1 + n2), sep = "")
      outcomes <- c(rep("A", n1), rep("B", n2))
      names(outcomes) <- colnames(data)
      results <-
        EMDomics::calculate_emd(
          data,
          outcomes,
          nperm = 10,
          parallel = FALSE,
          verbose = 0
        )
      for (i in 1:p) {
        w[k, i] <- results$emd[i, 2]
      }
    }
  }

  w <- colMins(w, value = TRUE)
  w <- p * exp(w) / sum(exp(w))
  return(w)
}

CDN_fast <-
  function(X,
           lambda1 = NULL,
           lambda2 = NULL,
           err_threshold = 1e-6,
           truncate = 1e-4,
           verbose = 1,
           ifparallel = FALSE) {
    library(EMDomics)
    library(Rfast)
    library("progress")
    if (ifparallel) {
      library("parallel")
      library("foreach")
      library("iterators")
      library("doParallel")
      cl <- makeCluster(20)
      registerDoParallel(cl)
    }
    # X: matrix list shape:(K,2),each one shape:(n_ki,p)
    K <- length(X) / 2
    p <- dim(X[[1]])[2]
    n <- matrix(unlist(lapply(X, function(x) {
      dim(x)[1]
    })), nrow = K)
    print("Intergrating p-values......")
    w <- calculate_w(X, K, p, ifparallel = ifparallel)
    if (is.null(lambda1)) {
      lambda1 <- 0.1 * mean(n)
    }
    if (is.null(lambda2)) {
      lambda2 <- 1 * mean(n)
    }

    sigma_k1 <-
      array(unlist(lapply(X[1:K], function(x) {
        cov(x)
      })), c(p, p, K))
    sigma_k2 <-
      array(unlist(lapply(X[(K + 1):(2 * K)], function(x) {
        cov(x)
      })), c(p, p, K))
    sigma_k1 <- aperm(sigma_k1, c(3, 1, 2))
    sigma_k2 <- aperm(sigma_k2, c(3, 1, 2))
    I <- diag(rep(1, p))
    I_k <- array(rep(I, K), c(p, p, K))
    I_k <- aperm(I_k, c(3, 1, 2))
    # Initialization primiary variables
    theta1 <- I_k
    theta2 <- I_k
    theta1_old <- I_k
    theta2_old <- I_k

    # Auxiliary Variables
    Z1 <- I_k
    Z2 <- I_k
    V <- I_k
    W <- I_k

    # Dual Variables
    O <- array(rep(0, K * p * p), c(K, p, p))
    F <- O
    G <- O
    Q1 <- O
    Q2 <- O


    # parameters
    rho0 <- 10 # larger for faster to coverge
    miu <- 2
    rho_max <- 5**3
    step_max <- 1000 # much enough steps to ensure convergency
    rho <- rho0

    # progress bar
    proceed.now <- 0
    pb <- progress_bar$new(
      format = "  Solving [:bar] :percent worst eta: :eta",
      total = (step_max),
      clear = FALSE,
      width = 80
    )
    print("Start Solving!")
    start <- Sys.time()
    for (step in 1:step_max) {
      pb$tick()
      # Primary Variables
      A1_ <- (n[, 1] * sigma_k1
        + F
        + Q1
        - rho * (theta2_old + V + W + Z1)) / n[, 1]
      theta1 <- update_theta_k(A1_, 2 * rho / n[, 1], ifparallel = ifparallel)

      A2_ <- (n[, 2] * sigma_k2
        - F
        + Q2
        + rho * (-theta1_old + V + W - Z2)) / n[, 2]
      theta2 <- update_theta_k(A2_, 2 * rho / n[, 2], ifparallel = ifparallel)

      # Auxiliary Variables
      Z1 <- soft_threshold(rho * theta1 + Q1, lambda1) / rho
      Z2 <- soft_threshold(rho * theta2 + Q2, lambda1) / rho
      C_ <- F + rho * (theta1 - theta2 - W + aperm(W, c(1, 3, 2)))
      V <- update_V_new(C_, rho, lambda2 * w)
      W <- 0.5 * (theta1 + aperm(V, c(1, 3, 2)) - V - theta2) + 0.5 / rho * (F + aperm(G, c(1, 3, 2)))

      # Dual Variables
      F <- F + rho * (theta1 - theta2 - (V + W))
      G <- G + rho * (V - aperm(W, c(1, 3, 2)))
      Q1 <- Q1 + rho * (theta1 - Z1)
      Q2 <- Q2 + rho * (theta2 - Z2)

      err <- max(
        sqrt(sum((
          theta1 - theta1_old
        )^2)) / sqrt(sum((theta1_old)^2)),
        sqrt(sum((
          theta2 - theta2_old
        )^2)) / sqrt(sum((theta2_old)^2))
      )

      if (err <= err_threshold) {
        break
      }

      if (step %% 100 == 0 & verbose == 1) {
        print(err)
      }
      # obj = objective(theta1,theta2,S1,S2,V,lambda1,lambda3,n)
      theta1_old <- theta1
      theta2_old <- theta2
      rho <- rho * miu
    }
    if (ifparallel) {
      stopCluster(cl)
    }
    end <- Sys.time()
    used_time <- end - start
    theta1 <- (theta1 + aperm(theta1, c(1, 3, 2))) / 2
    theta2 <- (theta2 + aperm(theta2, c(1, 3, 2))) / 2
    Z1 <- (Z1 + aperm(Z1, c(1, 3, 2))) / 2
    Z2 <- (Z2 + aperm(Z2, c(1, 3, 2))) / 2
    theta1[abs(theta1) < truncate] <- 0
    theta2[abs(theta2) < truncate] <- 0
    Z1[abs(Z1) < truncate] <- 0
    Z2[abs(Z2) < truncate] <- 0
    V[abs(V) < truncate] <- 0
    if (verbose == 1 & step == step_max & rho >= rho_max) {
      print(paste0(
        "Algorithm done with maximum iterations! relative error = ",
        round(err, 6)
      ))
    } else if (verbose == 1) {
      print(paste0("Algorithm has converged! relative error = ", round(err, 6)))
    }
    partial1 <- theta2partial(theta1)
    partial2 <- theta2partial(theta2)
    delta <- partial2 - partial1
    delta[abs(delta) < truncate] <- 0
    return(list(
      delta = delta,
      theta1 = theta1,
      theta2 = theta2,
      V = V,
      W = W,
      Z1 = Z1,
      Z2 = Z2
    ))
  }