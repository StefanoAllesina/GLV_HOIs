source("equilibria_simplified_model.R")
source("dynamics_simplified_model.R")
source("plotting.R")
set.seed(3) # interesting
set.seed(9) # stable middle
#set.seed(6) # bistable
#set.seed(19) # tristable (2spp)
pars <- rnorm(4)
# only competitive interactions
pars <- - abs(pars)
pars[1] <- abs(pars[1])
eq <- get_equilibria(pars)
dynamics <- integrate_multiple_starting_points(8, pars)
# plot biomass
show(pl_biomass(dynamics, pars, eq))
# plot V
show(pl_vstar(dynamics, pars, eq))
# plot dynamics
show(pl_dynamics_tern(dynamics, pars, eq))
print(eq)

