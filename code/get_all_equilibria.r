#Calculate all equilibria of a colection of GLV systems with HOIs

#local paths to be changed by user
setwd("/Users/pablolechon/Desktop/phd/GLV_HOIs/code")
#find the path to the julia executable in your machine by running the following 
#command in a julia REPL
#julia> println(Sys.BINDIR)
path_to_julia = "/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia"
#path_to_julia = "/Users/pablo/julia-1.8.3/bin/julia"

#load auxiliary functions to build parameter sets
source("build_parameter_sets.r")

#specify simulation number and diversity levels
simulations = 2
diversities = rep(c(3, 4), each = simulations)

#stack parameters sets of different communities
pars = stack_parameters(diversities)

#get subcommunities of each community
pars_subcomm = subcomm_expand(pars, diversities)

#save parameter sets
write.table(pars_subcomm[[1]], "../data/rs.csv", col.names = F, row.names = F)
write.table(pars_subcomm[[2]], "../data/As.csv", col.names = F, row.names = F)
write.table(pars_subcomm[[3]], "../data/Bs.csv", col.names = F, row.names = F)
write.table(pars_subcomm[[4]], "../data/diversities.csv", col.names = F, row.names = F)

#find all roots of the induced system of polynomials through homotopy continuation
system(paste(path_to_julia, "find_roots.jl"))

#create column names
spp_names = unlist(lapply(seq(max(diversities)), paste, 'spp', sep = ""))
col_names = c('solution_ind', 'diversity', spp_names)

#load all equilibria
all_roots = read.csv("../data/roots.csv", sep = '\t', header = F,
                     col.names = col_names, check.names = F)

#merge with rest of parameter data 
data = c(all_roots, data.frame('simulation' = ))