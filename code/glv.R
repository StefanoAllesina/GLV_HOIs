library(deSolve) # to integrate ODEs
library(tidyverse) # plotting
THRESH <- 10^-10 # consider extinct if below the threshold

# parameter structure:
# r: vector of growth rates (n)
# A: matrix of pairwise interactions (n x n)
# B: list of matrices for HOIs (n x (n x n))

glv <- function(t, x, pars){
  # zero extinct species
  x[x < THRESH] <- 0
  dx <- x * (pars$r + as.vector(pars$A %*% x))
  return(list(dx))
}

glv_hois <- function(t, x, pars){
  # zero extinct species
  x[x < THRESH] <- 0
  # build vector HOIs
  HOIs <- unlist(lapply(pars$B, function(Bi) x %*% Bi %*% x))
  dx <- x * (pars$r + as.vector(pars$A %*% x) + HOIs)
  return(list(dx))
}

generate_random_pars <- function(n){
  # build a GLV with pairwise interactions that is feasible
  # with random pairwise interactions
  x <- abs(rnorm(n))
  # choose a stable equilibrium
  success <- FALSE
  while(!success){
    A <- matrix(rnorm(n * n), n, n)
    # make diagonal negative
    diag(A) <- -abs(diag(A))
    eA <- eigen(diag(x) %*% A, only.values = TRUE, symmetric = FALSE)$values
    if (max(Re(eA)) < 0) success <- TRUE
  }
  r <- as.vector(-A %*% x)
  # now build the tensor such that the equilibrium is unchanged
  # use only "interaction modification" 
  # for three variables, the only option is no HOIs
  B <- list()
  xxt <- x %o% x
  for (i in 1:n){
    Bi <- matrix(rnorm(n * n), n, n)
    diag(Bi) <- 0
    Bi[i,] <- 0
    Bi[,i] <- 0
    # subtract mean
    Bi[Bi !=0] <- Bi[Bi !=0] - mean(Bi[Bi !=0])
    # divide each element
    Bi <- Bi / xxt
    B[[i]] <- Bi
  }
  return(pars = list(
    n = n,
    r = r,
    A = A,
    B = B,
    xstar = x
  ))
}

test_pars <- function(n){
  # build random model
  tmp <- generate_random_pars(n)
  # growth rates at equilibrium (pairwise model)
  print(round(unlist(glv(0, tmp$xstar, tmp)), 15))
  # growth rates at equilibrium (HOIs model)
  print(round(unlist(glv_hois(0, tmp$xstar, tmp), 15)))
}

integrate_dynamics <- function(pars, maxtime = 100, steptime = 0.1){
  time_integrate <- seq(0, maxtime, by = steptime)
  # initial conditions are a perturbation of the equilibrium
  x0 <- abs(pars$xstar * (1 + rnorm(pars$n, 0, 0.001)))
  out_pairs <- ode(y = x0, times = time_integrate, func = glv, 
                   parms = pars, method = "ode45")
  out_hois <- ode(y = x0, times = time_integrate, func = glv_hois, 
                   parms = pars, method = "ode45")
  return(list(pars = pars,
              x0 = x0,
              out_pairs = out_pairs,
              out_hois = out_hois))
}

plot_output <- function(out){
  dt <- out$out_pairs %>% 
    as.data.frame() %>% 
    pivot_longer(names_to = "variable", values_to = "abundance", cols = -time) %>% 
    add_column(model = "pairs")
  dt <- dt %>% bind_rows(out$out_hois %>% 
                           as.data.frame() %>% 
                           pivot_longer(names_to = "variable", values_to = "abundance", cols = -time) %>% 
                           add_column(model = "hois"))
    pl <- ggplot(dt, aes(x = time, y = abundance, colour = variable, linetype = model)) + 
      geom_line() + facet_wrap(~variable, scales = "free") + theme_bw()
    return(pl)
}

test_integration <- function(n){
  pars <- generate_random_pars(n)
  out <- integrate_dynamics(pars)
  show(plot_output(out))
  return(out)
}