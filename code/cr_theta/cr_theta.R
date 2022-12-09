## GOAL:
## This code shows that for the theta-logistic
## consumer-resource model, feasible equilibria
## are globally stable
## Call check_global(3, 0.5, 100)
## to simulate from 100 initial conditions
## the dynamics of the model when there are 
## 3 consumers and 3 resources and theta = 1/2
## the plot shows the (planted) feasible 
## equilibrium in red

library(deSolve)
library(tidyverse)
THRESH <- 10^-10

build_pars <- function(n, theta){
  # build a theta-logistic CR model
  # with a feasible equilibrium
  ys <- runif(n)
  xs <- runif(n)
  C <- matrix(runif(n * n), n, n)
  # resource growth rates
  r <- as.vector(ys^theta + t(C) %*% xs)
  # consumers death rates
  d <- as.vector(C %*% ys)
  return(list(
    r = r,
    d = d,
    C = C,
    theta = theta,
    xs = xs,
    ys = ys,
    n = n
  ))
}

cr_theta <- function(time, z, pars){
  z[z < THRESH] <- 0
  x <- z[1:pars$n]
  y <- z[pars$n + 1:pars$n]
  dx <- x * as.vector(-pars$d + pars$C %*% y)
  dy <- y * as.vector(pars$r - y^pars$theta - t(pars$C) %*% x)
  return(list(c(dx, dy)))
}

test_pars <- function(n, theta){
  pr <- build_pars(n, theta)
  print(pr)
  # check equil condition
  print(cr_theta(0, c(pr$xs, pr$ys), pr))
}

integrate_dynamics <- function(pars, z0 = NULL, maxint = 1000, stepint = 1){
  if (is.null(z0)) z0 <- runif(2 * c(pars$xs, pars$ys))
  out <- ode(y = z0, 
             times = seq(0, maxint, by = stepint), 
             func = cr_theta, parms = pars, method = "ode45")
  return(out)
}

# check global stability
check_global <- function(n, theta, tries){
  pr <- build_pars(n, theta)
  output <- tibble(spp = 1:(2*n),
                   value = c(pr$xs, pr$ys),
                   type = "equil",
                   rep = 0
                   )
  for (i in 1:tries){
    out <- integrate_dynamics(pr)
    output <- bind_rows(output,
                        tibble(
                          spp = 1:(2*n),
                          value =as.numeric(tail(out,1)[-1]),
                          type = "dynamics",
                          rep = i
                        ))
  }
  plgs <- ggplot(output %>% filter(type == "dynamics")) + 
    aes(x = as.factor(spp), y = value) + geom_boxplot() + 
    geom_point(data = output %>% filter(type == "equil"), colour = "red")
  return(plgs)
}
