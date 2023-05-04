#Calculate all equilibria of a colection of GLV systems with HOIs where parameters are sampled
#with different constraints

#find the path to the julia executable in your machine by running the following 
#command in a julia REPL
#julia> println(Sys.BINDIR)
path_to_julia = "/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia" #for my mac
#path_to_julia = "/home/pablo/julia-1.8.3/bin/julia" #for my desktop

#load auxiliary functions to build parameter sets
source("build_parameter_sets.r")

#specify simulation number and diversity levels
simulations = 1000
diversities = c(3, 4, 5, 6, 7, 8, 9, 10)
div_sim = rep(diversities, each = simulations)
par_type = 'random'

#sample and stack parameters sets of different communities
pars = sample_stack_par(div_sim, par_type)#A and B are symmetric
#get subcommunities of each community
#pars_subcomm = subcomm_expand(pars, div_sim)
#save parameter sets
save_parameter_set(pars, par_type)
#find all roots of the induced system of polynomials through homotopy continuation
system(paste(path_to_julia, "find_roots.jl", par_type))
#create column names and save file
spp_names = seq(max(diversities))
col_names = c('pool_div', 'pool_id', 'par_set_id',  spp_names)
all_roots = read.csv(paste("../data/roots_", par_type, ".csv", sep = ""), 
                     sep = '\t', header = F,
	       	     col.names = col_names, check.names = F)
#save clean data
write.csv(all_roots, paste("../data/roots_", par_type, ".csv", sep = ""), 
	  row.names = F)

