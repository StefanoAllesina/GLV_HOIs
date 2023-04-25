#Calculate all equilibria of a colection of GLV systems with HOIs where parameters are sampled
#with different constraints

#find the path to the julia executable in your machine by running the following 
#command in a julia REPL
#julia> println(Sys.BINDIR)
#path_to_julia = "/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia" #for my mac
path_to_julia = "/home/pablo/julia-1.8.3/bin/julia" #for my desktop

#load auxiliary functions to build parameter sets
source("build_parameter_sets.r")

#specify simulation number and diversity levels
simulations = 10
diversities = c(3, 6, 9)
div_sim = rep(diversities, each = simulations)
par_type = 'symmetric'

#sample and stack parameters sets of different communities
pars = sample_stack_par(div_sim, par_type)#A and B are symmetric
#get subcommunities of each community
pars_subcomm = subcomm_expand(pars, div_sim)
#save parameter sets
save_parameter_set(pars_subcomm, par_type)
#find all roots of the induced system of polynomials through homotopy continuation
system(paste(path_to_julia, "find_roots.jl", par_type))
#create column names
spp_names = unlist(lapply(seq(max(diversities)), paste, 'spp', sep = ""))
col_names = c('solution_id', 'diversity', 'pool',  spp_names)
#load all equilibria
all_roots = read.csv(paste("../data/roots_", par_type, ".csv", sep = ""), 
                     sep = '\t', header = F,
                     col.names = col_names, check.names = F)
