build_GLV_HOIs <- function(n, mode = "stable", HOIs = "modification", zerosumBi = TRUE){
  # step 1:
  # find a feasible GLV that is locally stable, neutrally stable, or unstable
  x <- abs(rnorm(n))
  success <- FALSE
  while(!success){
    if (mode == "stable"){
      A <- matrix(rnorm(n * n), n, n)
      # if trace is positive, make negative
      if (sum(diag(A)) > 0) diag(A) <- -diag(A)
      B <- diag(x) %*% A
      Rel1B <- max(Re(eigen(B, only.values = TRUE)$values))
      if (Rel1B < 0) success <- TRUE
    }
    if (mode == "neutral"){
      # take A to be skew-symmetric (only sensible for n even)
      A <- matrix(rnorm(n * n), n, n)
      A <- A - t(A)
      success <- TRUE
    }
    if (mode == "unstable"){
      A <- matrix(rnorm(n * n), n, n)
      B <- diag(x) %*% A
      Rel1B <- max(Re(eigen(B, only.values = TRUE)$values))
      if (Rel1B > 0) success <- TRUE
    }
  }
  r <- as.vector(-A %*% x)
  # step 2: 
  # add Higher-order interactions
  # case 1: adding HOIs to GLV does not change equilibrium
  if (zerosumBi){
    B <- list()
    xxt <- x %o% x
    for (i in 1:n){
      Bi <- matrix(rnorm(n * n, 0, 0.1), n, n)
      if (HOIs == "modification"){
        diag(Bi) <- 0
        Bi[i,] <- 0
        Bi[,i] <- 0
      }
      # make upper-triangular
      Bi <- Bi + t(Bi) - diag(diag(Bi))
      Bi[lower.tri(Bi)] <- 0
      # subtract mean
      Bi[Bi !=0] <- Bi[Bi !=0] - mean(Bi[Bi !=0])
      # divide each element
      Bi <- Bi / xxt
      B[[i]] <- Bi
    } 
  } else{
      # in this case, any weighted average of pairwise and hois does not modify equilibrium
      B <- list()
      xxt <- x %o% x
      for (i in 1:n){
        Bi <- matrix(rnorm(n * n, 0, 0.1), n, n)
        if (HOIs == "modification"){
          diag(Bi) <- 0
          Bi[i,] <- 0
          Bi[,i] <- 0
        }
        # make upper-triangular
        Bi <- Bi + t(Bi) - diag(diag(Bi))
        Bi[lower.tri(Bi)] <- 0
        # subtract mean
        Bi[Bi !=0] <- Bi[Bi !=0] - mean(Bi[Bi !=0])
        # add negative growth rates
        Bi[Bi !=0] <- Bi[Bi !=0] - r[i] / sum(Bi !=0)
        # divide each element
        Bi <- Bi / xxt
        B[[i]] <- Bi
    }
  } 
  return(pars = list(
    n = n,
    r = r,
    A = A,
    B = B,
    xstar = x
  ))
}

test_pars <- function(n=4){
  print("Test stable, modification, zerosumB")
  # build random model
  p1 <- build_GLV_HOIs(n = n, mode = "stable", HOIs = "modification", zerosumBi = TRUE)
  # growth rates at equilibrium (pairwise model)
  print(round(unlist(glv(0, p1$xstar, p1)), 15))
  # growth rates at equilibrium (HOIs model)
  print(round(unlist(glv_hois(0, p1$xstar, p1), 15)))
  
  print("Test neutrally stable, modification, zerosumB")
  # build random model
  p2 <- build_GLV_HOIs(n = n, mode = "neutral", HOIs = "modification", zerosumBi = TRUE)
  # growth rates at equilibrium (pairwise model)
  print(round(unlist(glv(0, p2$xstar, p2)), 15))
  # growth rates at equilibrium (HOIs model)
  print(round(unlist(glv_hois(0, p2$xstar, p2), 15)))
  
  print("Test unstable, modification, zerosumB")
  # build random model
  p3 <- build_GLV_HOIs(n = n, mode = "unstable", HOIs = "modification", zerosumBi = TRUE)
  # growth rates at equilibrium (pairwise model)
  print(round(unlist(glv(0, p3$xstar, p3)), 15))
  # growth rates at equilibrium (HOIs model)
  print(round(unlist(glv_hois(0, p3$xstar, p3), 15)))
  
  print("Test stable, modification, general B")
  # build random model
  p1 <- build_GLV_HOIs(n = n, mode = "stable", HOIs = "general", zerosumBi = TRUE)
  # growth rates at equilibrium (pairwise model)
  print(round(unlist(glv(0, p1$xstar, p1)), 15))
  # growth rates at equilibrium (HOIs model)
  print(round(unlist(glv_hois(0, p1$xstar, p1), 15)))
  
  print("Test neutrally stable, modification, general B")
  # build random model
  p2 <- build_GLV_HOIs(n = n, mode = "neutral", HOIs = "general", zerosumBi = TRUE)
  # growth rates at equilibrium (pairwise model)
  print(round(unlist(glv(0, p2$xstar, p2)), 15))
  # growth rates at equilibrium (HOIs model)
  print(round(unlist(glv_hois(0, p2$xstar, p2), 15)))
  
  print("Test unstable, modification, general B")
  # build random model
  p3 <- build_GLV_HOIs(n = n, mode = "unstable", HOIs = "general", zerosumBi = TRUE)
  # growth rates at equilibrium (pairwise model)
  print(round(unlist(glv(0, p3$xstar, p3)), 15))
  # growth rates at equilibrium (HOIs model)
  print(round(unlist(glv_hois(0, p3$xstar, p3), 15)))
  
  print("Test stable, modification, mix")
  # build random model
  p1 <- build_GLV_HOIs(n = n, mode = "stable", HOIs = "modification", zerosumBi = FALSE)
  # now produce different combinations
  p10 <- p1
  p10$B <- lapply(p1$B, function(m) m * 0)
  p10$A <- p1$A * 1
  print("Only pairs")
  print(round(unlist(glv_hois(0, p10$xstar, p10), 15)))
  p01 <- p1
  p01$B <- lapply(p1$B, function(m) m * 1)
  p01$A <- p1$A * 0
  print("Only HOIs")
  print(round(unlist(glv_hois(0, p01$xstar, p01), 15)))
  ph <- p1
  ph$B <- lapply(p1$B, function(m) m * 0.5)
  ph$A <- p1$A * 0.5
  print("50% pairs, 50% HOIs")
  # growth rates at equilibrium (HOIs model)
  print(round(unlist(glv_hois(0, ph$xstar, ph), 15)))
}

