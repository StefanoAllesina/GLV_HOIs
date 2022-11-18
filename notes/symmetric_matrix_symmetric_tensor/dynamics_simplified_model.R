library(deSolve)
THRESH <- 10^-10

sm <- function(t, x, pars){
  x[x < THRESH] <- 0
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  r <- pars$r
  a1 <- pars$a1
  a2 <- pars$a2
  b <- pars$b
  dx1 <- x1 * (r + a1 * x1 + a2 * x2 + a2 * x3 + b * x2 * x3)
  dx2 <- x2 * (r + a2 * x1 + a1 * x2 + a2 * x3 + b * x1 * x3)
  dx3 <- x3 * (r + a2 * x1 + a2 * x2 + a1 * x3 + b * x1 * x2)
  return(list(c(dx1, dx2, dx3)))
}

integrate_sm <- function(pars, x0 = NULL, maxtime = 100, bytime = 0.1){
  time_integrate <- seq(0, maxtime, by = bytime)
  parameters <- list(r = pars[1],
                     a1 = pars[2],
                     a2 = pars[3],
                     b = pars[4])
  # random initial conditions
  if (is.null(x0)) x0 <- abs(rnorm(3) * -2*pars[1]/pars[2])
  out <- ode(y = x0, times = time_integrate, func = sm, 
            parms = parameters, method = "ode45")
  return(list(pars = pars,
              x0 = x0,
              out = out))
}

integrate_multiple_starting_points <- function(nruns = 3, pars, maxtime = 100, bytime = 0.1){
  output <- tibble()
  for (i in 1:nruns){
    out <- integrate_sm(pars, NULL, maxtime, bytime)$out
    out <- out %>% 
      as.data.frame() %>% 
      pivot_longer(names_to = "variable", values_to = "abundance", cols = -time) %>% 
      add_column(run = letters[i])
    output <- bind_rows(output, out)
  }
  output <- output %>% mutate(abundance = ifelse(abundance < THRESH, 0, abundance))
  return(output)
}
