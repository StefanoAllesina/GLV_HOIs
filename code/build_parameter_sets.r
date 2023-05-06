#Set of functions used in get_all_equilibria.r to build parameter sets of
#different communities and, if desired, of all possible subcommunities 
#from each community. The goal is to build a data strucutre that can be 
#used by Julia to find all the roots at once. Additionally, there are also
#functions to sample tensor B with different constraints. All samples are
#done from a normal distribution with mean 0 and sigma 0.1

library(gtools)
source("generate_B.R") #sample symmetric and constrained B

all_presence_combs = function(n){
  #get all possible species combinations
  #n = 3 will give (1), (2), (3), (1,2), (1,3), (2,3), (1,2,3)
  
  v = seq(n)
  combs = do.call("c", lapply(seq_along(v), function(i) combn(v, i, FUN = list)))
  return(combs)
}

zeroed = function(r, A, B, spp_vec){
  #Zero elements in r, A, and B corresponding to absent species (spp_vec)
  
  n = length(which(r!=0))
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
          B[(i-1)*n + j,k] = 0 #transform three indices to two stacked indices
        }
      }
    }
  }
  return(list(r, A, B))
}

build_pars_subcomm = function(n, r, A, B){
  #given parameters of a species poll, build the parameter sets of 
  #all 2^n - 1 sub-communities
  
  #get all combinations of species
  spp_combs = all_presence_combs(n)
  n_sub = 2^n-1 #total number of subcommunities
  #place holder for subcommunity diversity
  n_spp_vec = rep(0, n_sub)
  #preallocate long matrices to hold all parameters combinations
  rs = matrix(0, nrow = n_sub, ncol = n)
  As = matrix(0, nrow = n_sub*n, ncol = n)
  Bs = matrix(0, nrow = n_sub*n^2, ncol = n)
  #loop over subcommunities
  for (i in 1:n_sub){
    #get ith species combination 
    ind = spp_combs[[i]]
    #zero out such indices in r, A, B
    pars = zeroed(r, A, B, ind)
    #parse
    r_sub = pars[[1]]
    A_sub = pars[[2]]
    #for B, transform list of lists to a long matrix prior to parsing
    B_sub = matrix(unlist(pars[[3]]), nrow = n^2, ncol = n)
    #number of species in subcommunity
    n_spp_vec[i] = length(ind)
    #allocate to holder long matrices
    rs[i,] = r_sub
    As[(n*(i-1)+1):(n*i),] = A_sub
    Bs[(n^2*(i-1)+1):(n^2*i),] = B_sub
  }
  return(list(rs, As, Bs, n_spp_vec))
}

sample_parameter_set = function(n, constraints = 'random',
				stability = 'random'){
  #sample parameter set with n species
  #kept as separate function to implement constraints modularly
  
  #sample sgrowth rates
  r = rnorm(n, 0, 0.1)
  #sample A and B from uniform distribution
  if (constraints == 'random'){
    #interaction matrix
    A = matrix(rnorm(n^2, 0, 0.1), n, n)
    #tensor of HOIs
    B = list()
    B[[1]] = matrix(rnorm(n^2, 0, 0.1), n, n)
    B[[2]] = matrix(rnorm(n^2, 0, 0.1), n, n)
    B[[3]] = matrix(rnorm(n^2, 0, 0.1), n, n)
    B_mat = matrix(unlist(B), nrow = n^2, ncol = n)
  }
  #sample symmetric A and B symmetric from uniform distribution
  else if (constraints == 'symmetric'){
    A_1 = matrix(rnorm(n^2, 0, 0.1), n, n)
    A = 1/2*(A_1 + t(A_1))
    B_mat = matrix(unlist(sample_symmetric_B(n)), nrow = n^2, ncol = n)
  }
  #sample parameters to preserve equilibrium. any weighted average of 
  #pairwise and hois does not modify equilibrium (implemented by stefano)
  else if (constraints == 'prerserving_sa'){
    #sample feasible equilibrium to preserve
    x <- abs(rnorm(n))
    #initialize B
    B <- list()
    #get xx^T
    xxt <- x %o% x
    for (i in 1:n){
      Bi <- matrix(rnorm(n * n, 0, 0.1), n, n)
      # make upper-triangular
      Bi <- Bi + t(Bi) - diag(diag(Bi))
      Bi[lower.tri(Bi)] <- 0
      # subtract mean
      Bi[Bi !=0] <- Bi[Bi !=0] - mean(Bi[Bi !=0])
      # divide each element
      Bi <- Bi / xxt
      B[[i]] <- Bi
    }
  }
  #same case implemented by pablo
  else if (constraints == 'prerserving_pla'){
    #sample r
    r = runif(n)
    #get column sum vector
    row_sums = -alpha *r
    B = sample_preserving_B(row_sums, n)
  }
  return(list(r, A, B_mat))
}

get_indices = function(index, unit_history){
    #get from, to row indices for each parameter group
    #get the actual index weighted by the unit hitory
    from = sum(unit_history[1:index-1]) + 1
    to = from + unit_history[index] - 1
    return(c(from, to))
}

sample_stack_par = function(diversities, sampling_constr){
    #Sample parameter sets (r, A, B) according to constraints and 
    #stack them into a common table

    #get largest community and number of communities 
    n_max = max(diversities)
    n_com = length(diversities)
    #build a species pool vector, common to each community
    pool_vec = seq(n_com)
    #place holder for all parameter sets
    r_stacked = matrix(0, nrow = sum(n_com), ncol = n_max)
    A_stacked = matrix(0, nrow = sum(diversities), ncol = n_max)
    B_stacked = matrix(0, nrow = sum(diversities^2), ncol = n_max)
    #loop through all communities
    for (i in seq(n_com)){
    	n = diversities[i]
	    #sample parameters
	    par_set = sample_parameter_set(n, sampling_constr)
	    #store
	    r_stacked[i, 1:n] = par_set[[1]]
	    ind_A = get_indices(i, diversities)
	    A_stacked[ind_A[1]:(ind_A[2]), 1:n] = par_set[[2]]
	    ind_B = get_indices(i, diversities^2)
	    B_stacked[ind_B[1]:(ind_B[2]), 1:n] = par_set[[3]]
	  }
    return(list(r_stacked, A_stacked, B_stacked, diversities, pool_vec))
}

subcomm_expand = function(all_parameters, diversities){
  #expand parameter sets for each community in the stacked parameter sets
  #to include parameters of all subcommunities.
  
  #parse parameters families rs, As, Bs
  rs = all_parameters[[1]]
  As = all_parameters[[2]]
  Bs = all_parameters[[3]]
  #get size of biggest community and number of communities
  n_max = max(diversities)
  n_com = length(diversities)
  #number of subcommunities in each community
  n_sub_vec = 2^diversities-1
  #build a diversity vector matching dimensions of data with subcommunities
  diversities_sub = rep(diversities, n_sub_vec)
  #build a species pool vector, common to each community and corresponding
  #subcommunities
  pool_vec = rep(seq(n_com), times = 2^(diversities)-1)
  #place holder for all expanded parameter families
  r_expanded = matrix(0, nrow = sum(n_sub_vec), ncol = n_max)
  A_expanded = matrix(0, nrow = sum(n_sub_vec*diversities), ncol = n_max)
  B_expanded = matrix(0, nrow = sum(n_sub_vec*diversities^2), ncol = n_max)
  #index counters for efficient storing
  curr_ind_r = 1
  curr_ind_A = 1
  curr_ind_B = 1
  #loop through diversity levels
  for (i in seq(n_com)){
    #number of species and sub-communities in ith community
    n = diversities[i]
    n_sub = n_sub_vec[i]
    #get parameter set of ith community
    r = rs[i,1:n]
    ind_A = get_indices(i, diversities)
    A = As[ind_A[1]:(ind_A[2]),1:n]
    ind_B = get_indices(i, diversities^2)
    B = Bs[ind_B[1]:(ind_B[2]),1:n]
    #build parameter sets of all subcommunities of the ith community
    pars_subcomm = build_pars_subcomm(n, r, A, B)
    r_expanded[curr_ind_r:(curr_ind_r + n_sub-1), 1:n] = pars_subcomm[[1]]
    A_expanded[curr_ind_A:(curr_ind_A + n_sub*n-1), 1:n] = pars_subcomm[[2]]
    B_expanded[curr_ind_B:(curr_ind_B + n_sub*n^2-1), 1:n] = pars_subcomm[[3]]
    #update indices
    curr_ind_r = curr_ind_r + n_sub
    curr_ind_A = curr_ind_A + n_sub*n
    curr_ind_B = curr_ind_B + n_sub*n^2
  }
  return(list(r_expanded, A_expanded, B_expanded, diversities_sub, pool_vec))
}

save_parameter_set = function(parameter_set, output_name){
    #get how many parameter lists
    n_pars = length(parameter_set)
    for (i in seq(n_pars)){
	write.table(parameter_set[[i]], 
		    paste("../data/pars_", as.character(i), "_", output_name, ".csv", 
			  sep = ""),
		    col.names = F, row.names = F)
    }
}

load_parameter_set = function(path, number_parameters, par_suffix){
  pars = list()
  for (i in seq(number_parameters)){
    pars[[i]] = as.matrix(read.table(paste("../data/pars_", 
                               as.character(i),
                               '_',
                               par_suffix, 
                               ".csv",
                               sep = "" ),
                         header = F, sep = ' '))
  }
  return(pars)
}
