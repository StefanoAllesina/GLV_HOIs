using HomotopyContinuation #to solve system of polynomials
using LinearAlgebra #to compute dot products
using DelimitedFiles #to load and save files

function poly_i(x, r, A, B, i)
    #given parameters, construct system of polynomial equations
    r[i] + dot(A[i,:], x) + dot(x, B[:, :, i]*x)
end

#number of species (common to all communities)
n_spp = parse(Int64, ARGS[1])
#number of parameter sets
n_par_sets = 2^n_spp - 1
#change working directory
cd("/Users/pablolechon/Desktop/phd/GLV_HOIs/code") #change this to your own path
#read data
rs = readdlm("../data/rs.csv")
As = readdlm("../data/As.csv")
Bs = readdlm("../data/Bs.csv")
#preallocate storing array for results
global all_eq_mat = range(0,n_spp)|>collect 
#calculate equilibria for each subcommunity
for i in 1:n_par_sets
    #index ith parameter set
    r = rs[i, :]
    A = As[(n_spp*(i-1)+1):(n_spp*i), :]
    B = Bs[(n_spp^2*(i-1)+1):(n_spp^2*i), :]
    #transform B into a 3D array
    B = reshape(B, (n_spp, n_spp, n_spp))
    #declare dynamic variables
    @var x[1:n_spp];
    #initialize system of ODEs as empty list
    equations = [];
    #construct set of dynamic equations
    for j in 1:n_spp
        eqn = poly_i(x, r, A, B, j);
        append!(equations, eqn);
    end
    #build system to find solutions through homotopy
    F = System(equations)
    #solve system
    result = solve(F)
    #get only the real solutions 
    real_sols = real_solutions(result)
    #get number of real solutions
    n_sols = size(real_sols, 1)
    #create matrix of zeros with n_sols rows
    sol_mat = zeros(n_sols, n_spp)
    for k in 1:n_sols
        #get indices of present species
        inds = findall(x->x!=0, r)
        #assign equilibria to these positions
        sol_mat[k,inds] = real_sols[k]
    end
    #create column vector with simulation number
    sim = repeat([i], n_sols)
    #create rows to append to final data frame
    append_rows = hcat(sim, sol_mat)
    #store
    global all_eq_mat = hcat(all_eq_mat, append_rows')
end
#save them (by rows)
writedlm("../data/roots.csv", all_eq_mat')