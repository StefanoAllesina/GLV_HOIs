jacobian_index = function(xstar, A, B){
  #Calculate jacobian using index notation (inefficient)
  n = length(xstar)
  J = matrix(0, n, n)
  for (i in 1:n){
    for (j in 1:n){
      J[i, j] = A[i,j]*xstar[i] + xstar[i]*sum((B[[i]][j,]+B[[i]][,j])*xstar)
    }
  }
  return(J)
}

jacobian_vec = function(xstar, A, B){
  #######################################################
  #Calculate jacobian in a vectorized fashion (efficient)
  #Parameters:
    #xstar: Community equilibrium
    #A: Matrix of species interactions
    #B: Higher order interaction tensor
  #######################################################
  #linear part
  L = diag(xstar)%*%A
  #higher order interactions part
  H = sapply(1:n, function(i)(xstar[i]* ((B[[i]] + t(B[[i]])) %*% xstar) ))
  return(L + t(H))
}

################################################################################
#Timing each function
n_vec = seq(10, 500, 100)
n_n = length(n_vec)
#initialize time vectors
t_ind = rep(0, n_n)
t_vect = rep(0, n_n)
#calculate jacobian for increasing communnity size
for (i in seq(n_n)){
  n = n_vec[i]
  #sample params
  for (k in 1:n){
    B[[k]] = matrix(runif(n^2, 0,1), n, n)
  }
  xstar = runif(n, 0, 1)
  A = matrix(runif(n^2, 0, 1), n, n)
  #time functions
  t_0 = Sys.time()
  jacobian_index(xstar, A, B)
  t_f = Sys.time()
  t_ind[i] = as.numeric(t_f - t_0)
  t_0 = Sys.time()
  jacobian_vec(xstar, A, B)
  t_f = Sys.time()
  t_vect[i] = as.numeric(t_f-t_0)
}

#plot jacobian computation time vs community size for each function
plot(n_vec, t_ind, pch = 16)
points(n_vec, t_vect, pch = 16, col = 'blue')

