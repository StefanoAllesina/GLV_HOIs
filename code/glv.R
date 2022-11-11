library(deSolve) # to integrate ODEs
library(tidyverse) # plotting
source("generate_pars.R")
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

integrate_dynamics <- function(pars, model = "glv", alpha = 0, maxtime = 100, steptime = 0.1){
  time_integrate <- seq(0, maxtime, by = steptime)
  # initial conditions are a slight perturbation of the equilibrium
  x0 <- abs(pars$xstar * (1 + rnorm(pars$n, 0, 0.01)))
  if (model == "glv"){
    out <- ode(y = x0, times = time_integrate, func = glv, 
                   parms = pars, method = "ode45")
  } 
  if (model == "glv_hois"){
    out <- ode(y = x0, times = time_integrate, func = glv_hois, 
                     parms = pars, method = "ode45")
  }
  if (model == "mix"){
    pr2 <- pars
    pr2$A <- alpha * pars$A
    pr2$B <- (1 - alpha) * pars$B
    if (model == "glv_hois"){
      out <- ode(y = x0, times = time_integrate, func = glv_hois, 
                      parms = pars, method = "ode45")
    } 
  }
  return(list(pars = pars,
              x0 = x0,
              out = out))
}

plot_output <- function(out1, out2 = NULL){
  dt <- out1 %>% 
    as.data.frame() %>% 
    pivot_longer(names_to = "variable", values_to = "abundance", cols = -time) %>% 
    add_column(model = "first")
  if (!is.null(out2)){
    dt <- dt %>% bind_rows(out2 %>% 
                           as.data.frame() %>% 
                           pivot_longer(names_to = "variable", values_to = "abundance", cols = -time) %>% 
                           add_column(model = "second"))
  }
  pl <- ggplot(dt, aes(x = time, y = abundance, colour = variable, linetype = model)) + 
      geom_line() + facet_wrap(~variable, scales = "free") + theme_bw()
  return(pl)
}

test_integration <- function(n){
  pars <- build_GLV_HOIs(n = n, mode = "stable", HOIs = "modification", zerosumBi = TRUE)
  outpairs <- integrate_dynamics(pars, model = "glv")
  outhois <- integrate_dynamics(pars, model = "glv_hois")
  show(plot_output(outpairs$out, outhois$out))
  return(list(outpairs, outhois))
}