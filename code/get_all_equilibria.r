#Calculate all equilibria of a colection of GLV systems with HOIs

#local paths to be changed by user
setwd("/Users/pablolechon/Desktop/phd/GLV_HOIs/code")
path_to_julia = "/Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia"

#load auxiliary files
source("stack_params.r")
source("build_par_subcomms.r") #builds parameter sets of all sub-communities

#diversity levels and max level
n_vec = c(2,3)
n_sim = 2
#save parameter sets
write.table(rall, "../data/rs.csv", col.names = F, row.names = F)
write.table(Aall, "../data/As.csv", col.names = F, row.names = F)
write.table(Ball, "../data/Bs.csv", col.names = F, row.names = F)
df = data.frame(community, simulation, diversity)
write.table(df, "../data/parser.csv", col.names = F, row.names = F)
#find all roots of the induced system of polynomials through homotopy continuation
system(paste(path_to_julia, "find_roots.jl", as.character(n)))
#load all equilibria
all_roots = read.csv("../data/roots.csv", sep = '\t', header = T)
