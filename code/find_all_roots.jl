using HomotopyContinuation #to solve system of polynomials
using LinearAlgebra #to compute dot products
using DelimitedFiles #to load and save files

function poly_i(x, r, A, B, i)
    #given parameters, construct system of polynomial equations
    r[i] + dot(A[i,:], x) + dot(x, B[:, :, i]*x)
end

#change working directory
cd("/Users/pablolechon/Desktop/phd/GLV_HOIs/code")#change this to your path
#read data
r = readdlm("../data/r.csv")
A = readdlm("../data/A.csv")
B = readdlm("../data/B.csv")
#transform into a 3D array
n_spp = length(r)
B = reshape(B, (n_spp, n_spp, n_spp))

#declare variables
@var x[1:n_spp];
#initialize system as empty list
equations = [];
#construct collection with equations
for i in 1:n_spp
    eqn = poly_i(x, r, A, B, i);
    append!(equations, eqn);
end
#build system
F = System(equations)
#solve system
result = solve(F)
#get only the real solutions
real_sols = real_solutions(result)

#save them (by rows)
writedlm("../data/real_solutions.csv", real_sols)
