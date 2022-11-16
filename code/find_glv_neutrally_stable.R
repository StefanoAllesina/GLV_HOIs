source("glv.R")
# goal: find A, r such that 
# x* = -A^(-1) r > 0
# and D(x*)A has only imaginary eigenvalues
# one option A is skew symmetric, and has an even number of spp (to make it invertible)
n=3
test_integration_neutral <- function(n){
  pars <- build_GLV_HOIs(n, mode = "stable", HOIs = "else", zerosumBi = FALSE)
  sum_vec = -pars$r
  #get B 
  B = sample_B(sum_vec, n)
  pars$B = B
  outpairs <- integrate_dynamics(pars, model = "glv")
  outhois <- integrate_dynamics(pars, model = "glv_hois")
  show(plot_output(outpairs, outhois))
  return(out)
}

