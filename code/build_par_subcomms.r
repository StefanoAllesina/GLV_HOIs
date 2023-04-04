all_presence_combs = function(n){
  #get all possible combinations of present species
  v = seq(n)
  combs = do.call("c", lapply(seq_along(v), function(i) combn(v, i, FUN = list)))
  return(combs)
}

zeroed = function(r, A, B, spp_vec){
  #Zero r, A, and B elements of absent species
  for (i in 1:n){
    for (j in 1:n){
      for (k in 1:n){
        if (!(k %in% spp_vec)){
          r[k] = 0
        }
        if (!(j %in% spp_vec) | !(k %in% spp_vec)){
          A[j,k] = 0
        }
        if (!(i %in% spp_vec) | !(j %in% spp_vec) | !(k %in% spp_vec)){
          B[[i]][j,k] = 0
        }
      }
    }
  }
  return(list(r, A, B))
}

build_pars_subcomm(n, r, A, B){
  #get all combinations of species presence vectos
  present_combs = all_presence_combs(n)
  n_sub = 2^n-1
  #preallocate object that will hold all parameters
  rs = matrix(0, nrow = n_sub, ncol = n)
  As = matrix(0, nrow = n_sub*n, ncol = n)
  Bs = matrix(0, nrow = n_sub*n^2, ncol = n)
  #generate all parameter sets for all subcommunities
  for (i in 1:n_sub){
    #get ith vector of presence-absence
    ind = present_combs[[i]]
    pars = zeroed(r, A, B, ind)
    r_sub = pars[[1]]
    A_sub = pars[[2]]
    #transform list of lists to a long matrix
    B_sub = matrix(unlist(pars[[3]]), nrow = n^2, ncol = n)
    #assign to the bigger data structures
    rs[i,] = r_sub
    As[(n*(i-1)+1):(n*i),] = A_sub
    Bs[(n^2*(i-1)+1):(n^2*i),] = B_sub
  }
  return(list(rs, As, Bs))
}