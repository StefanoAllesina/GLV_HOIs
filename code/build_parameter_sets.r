#Set of functions used in get_all_equilibria.r to build parameter sets of
#different communities and, if desired, of all possible subcommunities 
#from each community. The goal is to build a data strucutre that can be 
#used by Julia to find all the roots at once.

#GET RID OF DIVERSITY CALCULATION IF NOT NEEDED
####################################################################
#maybe use this chunk later for merging roots and community parameters
####################################################################
##create indexing vectors to parse parameter sets efficiently in julia
#community = rep(rep(seq(n_com), n_sub_vec), each = n_sim)
#diversity = rep(rep(n_vec, n_sub_vec), each = n_sim)
#simulation = c()
#for (i in seq(n_com)){
#    vec = rep(seq(n_sim), each = n_sub_vec[i])
#    simulation = c(simulation, vec)
#}
##################################################################
#diversity junk code
##################################################################
#add diversity of subcommunities
#diversity[curr_ind_r:(curr_ind_r + n_sub - 1)] = pars_subcomm[[4]]

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

sample_parameter_set = function(n){
    #sample parameter set with n species
    #kept as separate function to implement constraints modularly

    #growth rates
    r = runif(n)
    #interaction matrix
    A = matrix(runif(n^2), n, n)
    #tensor of HOIs
    B = list()
    B[[1]] = -matrix(runif(n^2), n, n)
    B[[2]] = -matrix(runif(n^2), n, n)
    B[[3]] = -matrix(runif(n^2), n, n)
    B_mat = matrix(unlist(B), nrow = n^2, ncol = n)
    return(list(r, A, B_mat))
}

get_indices = function(index, unit_history){
    #get from, to row indices for each parameter group
    #get the actual index weighted by the unit hitory
    from = sum(unit_history[1:index-1]) + 1
    to = from + unit_history[index] - 1
    return(c(from, to))
}

stack_parameters = function(diversities){
    #Given a bunch of parameter sets (r, A, B), stack them into a common table

    #get largest community and number of communities 
    n_max = max(diversities)
    n_com = length(diversities)
    #place holder for all parameter sets
    r_stacked = matrix(0, nrow = sum(n_com), ncol = n_max)
    A_stacked = matrix(0, nrow = sum(diversities), ncol = n_max)
    B_stacked = matrix(0, nrow = sum(diversities^2), ncol = n_max)
    #loop through all communities
    for (i in seq(n_com)){
    	n = diversities[i]
	    #sample parameters
	    par_set = sample_parameter_set(n)
	    #store
	    r_stacked[i, 1:n] = par_set[[1]]
      ind_A = get_indices(i, diversities)
	    A_stacked[ind_A[1]:(ind_A[2]), 1:n] = par_set[[2]]
      ind_B = get_indices(i, diversities^2)
	    B_stacked[ind_B[1]:(ind_B[2]), 1:n] = par_set[[3]]
	  }
    return(list(r_stacked, A_stacked, B_stacked))
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
  return(list(r_expanded, A_expanded, B_expanded))
}