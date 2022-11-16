get_Q = function(sum_vec, n, tol = 1e-10, max_it = 100){
  ##############################################################################
  #Generate a symmetric matrix whose rows (columns) sum to sum_vec
  #Parameters
    #sum_vec (1xn array): array containing the row (column) sums
    #n (int): dimmension of the matrix
  #Returns
    #Q (nxn array)
  ##############################################################################
  #sample i.i.d. matrix
  M = matrix(runif(n^2), n, n)
  err = 1
  it = 0
  #make M doubly stochastic
  while (err>tol){
    row_sum = rowSums(M)
    #divide by row the row sums
    M_row = M/row_sum*sum_vec
    col_sum = colSums(M_row)
    #divide column by the row sums
    M = t(M_row)/col_sum*sum_vec
    #calculate error from a perfectly doubly stochastic matrix
    err = sum(abs(sum_vec-col_sum)) + sum(abs(sum_vec-row_sum))
    it = it + 1
    if (it > max_it){
      #non convergence, start over
      M = matrix(runif(n^2), n, n)
      it = 0
      err = 1
    }
  }
  Q = 1/2*(M+t(M))
  return(Q)
}
 s
get_B = function(Q, n){
  ##############################################################################
  #Given a matrix Q, compute a tensor B such that the row sums of each of its 
  #matrices sum to the elements of Q
  #Parameters:
    #Q (nxm array)
    #n (int)
  #Returns:
    #B (nxnxn list)
  ##############################################################################
  B = vector(mode = "list", length = n)
  for (i in seq(n)){
    for (j in seq(n)){
      Mi = matrix(runif(n^2), n,n)
      B[[i]] = Mi/rowSums(Mi)*Q[,i]
    }
  }
  return(B)
}

sample_B = function(row_sums, n){
  #get matrix Q
  Q = get_Q(row_sums, n)
  #get tensor B
  B = get_B(Q, n)
  return(B)
}

test = function(){
  #parameters
  n=4
  alpha = 0.5
  #sample r
  r = runif(n)
  #get column sum vector
  row_sums = -(1-alpha)*r
  B = sample_B(row_sums, n)
  return(B)
}

test()
