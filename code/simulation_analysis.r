#Script to analyze simulation results

#what simulations to analyze
par_type = 'symmetric'
#load all equilibria
all_roots = read.csv(paste("../data/roots_", par_type, ".csv", sep = ""), 
		     sep = '\t', header = F,
                     col.names = col_names, check.names = F)
#create a vector for subcommunities
