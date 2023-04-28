library(gtools)

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
  M = matrix(rnorm(n^2, 0, 0.1), n, n)
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
      M = matrix(rnorm(n^2, 0, 0.1), n, n)
      it = 0
      err = 1
    }
  }
  Q = 1/2*(M+t(M))
  return(Q)
}

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
      Mi = matrix(rnorm(n^2, 0, 0.1), n,n)
      B[[i]] = Mi/rowSums(Mi)*Q[,i]
    }
  }
  return(B)
}

sample_preserving_B = function(row_sums, n){
  #get matrix Q
  Q = get_Q(row_sums, n)
  #get tensor B
  B = get_B(Q, n)
  return(B)
}

is_permutation = function(vec1, vec2){
  #Check if vec1 is a permutation of vec2
  return(all(sort(vec1) == sort(vec2)))
}

perm2comb = function(permutation, combinations){
  #given a list of combinations and a permutation of one of these, identify 
  #which combination does permutation correspond to
  n_comb = nrow(combinations)
  for (i in seq(n_comb)){
    current_comb = combinations[i,]
    perm_is_comb = is_permutation(current_comb, permutation)
    if (perm_is_comb){
      return(i)
    }
  }
}

sample_symmetric_B = function(n){
  ##############################################################################
  #Build a random symmetric tensor
  #Parameters:
    #n (int)
  #Outputs:
    #B (nxnxn array)
  ##############################################################################
  #initialize 3D array for B
  B = array(rep(0, n^3), dim = c(n,n,n))
  #Get all the possible different indices combinations
  vec_comb = combinations(n, 3, repeats.allowed=TRUE)
  n_comb = nrow(vec_comb)
  #Assign random values to unique elements
  Bijk = rnorm(n_comb, 0, 0.1)
  for (i in seq(n)){
    for (j in seq(n)){
      for (k in seq(n)){
        #form a vector of indices
        current_permutation = c(i, j, k)
        #identify which of the combinations is this a permutation of
        ind_comb = perm2comb(current_permutation, vec_comb)
        #assign the element corresponding to that combination to this 
        #particular permutation
        B[i,j,k] = Bijk[ind_comb]
      }
    }
  }
  return(B)
  }

test = function(){
  #parameters
  n=4
  alpha = 0.5
  #sample r
  r = runif(n)
  #get column sum vector
  row_sums = -alpha *r
  B = sample_preserving_B(row_sums, n)
  B_symm = sample_symmetric_B(n)
  return(list(B, B_symm))
}

test()
