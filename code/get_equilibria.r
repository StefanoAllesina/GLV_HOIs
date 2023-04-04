setwd("/Users/pablolechon/Desktop/phd/GLV_HOIs/code")
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
#build parameter sets for all subcommunities
pars_subcomm = build_pars_subcomm(n, r, A, B)
#save
write.table(pars_subcomm[[1]], "../data/rs.csv", col.names = F, row.names = F)
write.table(pars_subcomm[[2]], "../data/As.csv", col.names = F, row.names = F)
write.table(pars_subcomm[[3]], "../data/Bs.csv", col.names = F, row.names = F)
#find all roots of the induced system of polynomials through homotopy continuation
system("/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia find_roots.jl 3")
#load all roots
all_roots = read.csv("../data/all_roots.csv", header = T)
