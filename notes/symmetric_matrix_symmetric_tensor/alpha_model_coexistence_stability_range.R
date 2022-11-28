source("equilibria_simplified_alpha_model.R")
source("dynamics_simplified_alpha_model.R")

set.seed(15)
pars <- rnorm(4)
# only competitive interactions
pars <- - abs(pars)
pars[1] <- abs(pars[1])


# To determine the range of alpha for which the coexistence equilibrium is stable
alpha <- seq(0, 1.0, 0.01)
res <- matrix(0,0,4) #stores the range of alpha for which the coexistence equilibria is stable

for (p in alpha){
  E <- get_equilibria(pars, p)
  for (i in E$label){
    if (all(sapply(E[i,1:3], function(x) x == E[i,1])) == TRUE & E[i,4] == TRUE){
      res <- rbind(res, c(as.numeric(E[i, 1:3]), p))
    }
  }
}
print(res)
