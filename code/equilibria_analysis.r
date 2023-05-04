#This file analyzes all equilibria of GLV HOIs system

library("tidyverse")
source("jacobian.R")
source("build_parameter_sets.r")

#set type of simulations to analyze
par_type = 'random'
#load data
data = read.csv(paste("../data/roots_", par_type, ".csv", sep = "" ))
rs = read.csv(paste("../data/pars_1_", par_type, ".csv", sep = ""), header = F, sep = ' ')
As = read.csv(paste("../data/pars_2_", par_type, ".csv", sep = ""), header = F, sep = ' ')
Bs = read.csv(paste("../data/pars_3_", par_type, ".csv", sep = ""), header = F, sep = ' ')
diversities = read.csv(paste("../data/pars_4_", par_type, ".csv", sep = ""), 
                       header = F, sep = ' ')
#get only indices of parameter sets that yielded at least one equilibrium
realized_comms = unique(data$par_set_id)
realized_diversities = diversities[realized_comms,]
#get dimmensions of data
n_cols = ncol(data)

#add useful metrics
df = data %>% add_count(par_set_id, name = 'n_eq') %>% #counts number of sols
  group_by(pool_id, pool_div, par_set_id) %>% 
  mutate(eq_id = seq(n_eq)) %>% #adds index to identify each solution
  select(-n_eq) %>% #delete unneded rows
  pivot_longer(cols = 4:n_cols, names_to = 'spp_id', 
               values_to = 'abundance') %>% 
  group_by(par_set_id, eq_id) %>% 
  mutate(comm_div = sum(abundance != 0), #gets sub-community abundance
         feasible = all(abundance >= 0)) #gets feasibility of equilibrium
  
#calculate stability of each equilibria
stability = function(eq, A, B){
  #calculate the jacobian
  J = jacobian_vec(eq, A, B)
  #get largest eigenvalue
  lambda_max = max(eigen(J, only.values = T))
  return(lambda_max)
}

matrix2list = function(matrix, dimmension){
  lst = lapply(split(seq_len(nrow(matrix)),(seq_len(nrow(matrix))-1) %/%dimension +1),
                        function(i) matrix[i,])
  return(lst)
}

all_stabilities = function(data){
  #preallocate vector of largest eigenvalue
  n_eq = nrow(data)
  n_col = ncol(data)
  lambda_max_vec = rep(0, )
  #separate data into groups, one for each subcommunity
  data_sep = data %>% group_by(par_set_id) %>% 
    group_split()
  n_subcomms = length(data_sep)
  #index for eigenvalue vector 
  idx = 1
  #loop over number of subcommunities
  for (i in seq(n_subcomms)){
    data_i = data_sep[[i]]
    n_sols = nrow(data_i)
    #loop over multiple equilibria
    for (j in seq(n_sols)){
      data_ij = data_i %>% filter(eq_id == j)
      #extract equilibrium
      eq_ij = data_ij$abundance
      #extract parameters
      par_set_ind = unique(data_ij$par_set_id)
      diversity_i = unique(data_ij$pool_div)
      ind_realized = which(realized_comms == par_set_ind)
      ind_A = get_indices(par_set_ind, realized_diversities)
      ind_B = get_indices(par_set_ind, realized_diversities^2)
      A_k = As[ind_A[1]:ind_A[2], 1:diversity_i]
      B_k = Bs[ind_B[1]:ind_B[2], 1:diversity_i]
      #transform B to list of lists to feed to jacobian function
      B = matrix2list(B_k, 3)
      lambda_max_vec[idx] = stability(eq_ij, A_k, B)
      idx = idx + 1
    }
  }
    
}
