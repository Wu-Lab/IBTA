library('parallel')
library('foreach')
library('iterators')
library('doParallel')
library('progress')

APC2 <- function(raw_X) {
  # X: (n, 2)
  # APC2: (n, 1)
  
  # protect from all same value
  raw_X[1, ] <- raw_X[1, ] + c(1e-12, 1e-12)
  tmp <- prcomp(raw_X, scale. = TRUE)
  APC2 <- abs(tmp$x[, 2])
  return(APC2)
}

pivot_APC2 <- function(raw_X,
                       pivot = F,
                       labels = NULL) {
  # X: (n, 2)
  # APC2: (n, 1)
  # labels: (n,1) ~ {0,1}
  ###### use pivot of median APC2
  
  
  if (!pivot) {
    raw_X[1, ] <- raw_X[1, ] + c(1e-12, 1e-12)
    tmp <- prcomp(raw_X, scale. = TRUE)
    APC2 <- abs(tmp$x[, 2])
    x_new <- NULL
    print("Calculate APC2 without pivot...")
    
  } else{
    if (is.null(labels)) {
      print("Need labels! Exit!")
      return()
    } else{
      raw_X0 <- raw_X[labels == 0, ]
      raw_X1 <- raw_X[labels == 1, ]
      
      raw_X0[1, ] <- raw_X0[1, ] + c(1e-12, 1e-12)
      tmp0 <- prcomp(raw_X0)
      APC2_0 <- abs(tmp0$x[, 2])
      
      raw_X1[1, ] <- raw_X1[1, ] + c(1e-12, 1e-12)
      tmp1 <- prcomp(raw_X1)
      APC2_1 <- abs(tmp1$x[, 2])
      
      if (median(APC2_0) > median(APC2_1)) {
        pivot <- 1
        tmp_pivot <- tmp1
        APC2_pivot <- APC2_1
        raw_X_new <- raw_X0
      } else{
        pivot <- 0
        tmp_pivot <- tmp0
        APC2_pivot <- APC2_0
        raw_X_new <- raw_X1
      }
      
      
      PC2_rotation <-
        scale(raw_X_new, center = tmp_pivot$center, scale = FALSE)
      PC2_rotation <- PC2_rotation %*% (tmp_pivot$rotation)
      APC2_new <- abs(PC2_rotation[, 2])
      APC2 <- raw_X[, 1]
      APC2[labels == pivot] <- APC2_pivot
      APC2[labels != pivot] <- APC2_new
      
      x_new <- raw_X
      x_new[labels == pivot, ] <- tmp_pivot$x
      x_new[labels != pivot, ] <- PC2_rotation
    }
  }
  return(list(APC2 = APC2, x_new = x_new))
}
