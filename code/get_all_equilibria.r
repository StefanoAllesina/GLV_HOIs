#Calculate all equilibria of a colection of GLV systems with HOIs

#local paths to be changed by user
setwd("~/Desktop/GLV_HOIs/code")
#find the path to the julia executable in your machine by running the following 
#command in a julia REPL
#julia> println(Sys.BINDIR)
path_to_julia = "/home/pablo/julia-1.8.3/bin/julia"

#load auxiliary functions to build parameter sets
source("build_parameter_sets.r")

#specify simulation number and diversity levels
simulations = 2
diversities = rep(c(3), each = simulations)

#build parameters for different communities
pars = stack_parameters(diversities)
diversities = rep(rep(diversities, n_sub_vec), each = n_sim)

#save parameter sets
write.table(pars[[1]], "../data/rs.csv", col.names = F, row.names = F)
write.table(pars[[2]], "../data/As.csv", col.names = F, row.names = F)
write.table(pars[[3]], "../data/Bs.csv", col.names = F, row.names = F)
write.table(diversities, "../data/diversities.csv", col.names = F, row.names = F)

#find all roots of the induced system of polynomials through homotopy continuation
system(paste(path_to_julia, "find_roots.jl", as.character(n)))

#load all equilibria
all_roots = read.csv("../data/roots.csv", sep = '\t', header = T)

#create column names
spp_names = unlist(lapply(seq(n), paste, 'spp', sep = ""))
col_names = c('sub_comm', spp_names)

#load all equilibria
all_roots = read.csv("../data/roots.csv", sep = '\t', header = T,
                     col.names = col_names, check.names = F)
