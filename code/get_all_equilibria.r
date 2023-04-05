#Example code on how to calculate the interior and all boundary equilibria of a GLV system with HOIs

#local paths to be changed by user
setwd("~/Desktop/GLV_HOIs/code")
#find the path to the julia executable in your machine by running the following 
#command in a julia REPL
#julia> println(Sys.BINDIR)
path_to_julia = "/home/pablo/julia-1.8.3/bin/julia"

#load auxiliary files
source("build_par_subcomms.r")

#number of species
n = 3
#growth rates
r = runif(n)
#interaction matrix
A = matrix(runif(n^2), n, n)
#tensor of HOIs
B = list()
B[[1]] = -matrix(runif(n^2), n, n)
B[[2]] = -matrix(runif(n^2), n, n)
B[[3]] = -matrix(runif(n^2), n, n)
#build ans save parameter sets for all subcommunities
build_pars_subcomm(n, r, A, B)
#find all roots of the induced system of polynomials through homotopy continuation
system(paste(path_to_julia, "find_roots.jl", as.character(n)))
#create column names
spp_names = unlist(lapply(seq(n), paste, 'spp', sep = ""))
col_names = c('sub_comm', spp_names)
#load all equilibria
all_roots = read.csv("../data/roots.csv", sep = '\t', header = T,
                     col.names = col_names, check.names = F)
