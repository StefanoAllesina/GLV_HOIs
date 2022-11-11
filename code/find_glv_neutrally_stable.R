source("glv.R")
# goal: find A, r such that 
# x* = -A^(-1) r > 0
# and D(x*)A has only imaginary eigenvalues
# one option A is skew symmetric, and has an even number of spp (to make it invertible)

test_integration_neutral <- function(n){
  pars <- generate_pars_neutral(n)
  out <- integrate_dynamics(pars)
  show(plot_output(out))
  return(out)
}

