source("equilibria_simplified_model.R")
source("dynamics_simplified_model.R")
source("plotting.R")


success <- FALSE
while(!success){
  pars <- rnorm(4)
  eq <- get_equilibria(pars)
  eq <- eq %>% mutate(nsp = (p1 > 0) + (p2 > 0) + (p3 > 0))
  nstable_coex <- sum((eq %>% filter(nsp == 3))$stability)
  if (nstable_coex > 1) {
    success <- TRUE
    print(eq)
  }
}