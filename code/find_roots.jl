using HomotopyContinuation #to solve system of polynomials
using LinearAlgebra #to compute dot products
using DelimitedFiles #to load and save files

function poly_i(x, r, A, B, i)
    #given parameters, construct glv_hoi system of polynomial equations
    r[i] + dot(A[i,:], x) + dot(x, B[:, :, i]*x)
end

#change working directory
cd("/Users/pablolechon/Desktop/phd/GLV_HOIs/code") #change this to your own path

#read data
#community index, simulation index, (sub)community diversity
parser = readdlm("../data/parser.csv")
#total number of communities
n_comms = size(parser,1)
#load growth rates, pairwise, three-way interactions
rs = readdlm("../data/rs.csv")
As = readdlm("../data/As.csv")
Bs = readdlm("../data/Bs.csv")
#place holder for storing all equilibria by rows
global all_eq_mat = range(0,max(n_spp_vec))|>collect 
#iterate through parameter sets (either communities or subcommunities)
for i in 1:n_comms
    #number of species in the ith (sub)community
    n_spp_i = parser[i,2]
    #get parameter set of ith (sub)community
    r = rs[i, :]
    A = As[(n_spp_i*(i-1)+1):(n_spp_i*i), :]
    B = Bs[(n_spp_i^2*(i-1)+1):(n_spp_i^2*i), :]
    #transform B into a 3D array
    B = reshape(B, (n_spp_i, n_spp_i, n_spp_i))
    #declare dynamic variables
    @var x[1:n_spp_i];
    #initialize system of ODEs as empty list
    equations = [];
    #construct set of dynamic equations
    for j in 1:n_spp_i
        eqn = poly_i(x, r, A, B, j);
        append!(equations, eqn);
    end
    #solve system and keep real solutions
    F = System(equations)
    result = solve(F)
    real_sols = real_solutions(result)
    #get number of solutions
    n_sols = size(real_sols, 1)
    #preallocate n_sols rows of zeros to store solutions
    sol_mat = zeros(n_sols, n_spp_i)
    #loop through solutions
    for k in 1:n_sols
        #use growth rate vector to find present species 
        inds = findall(x->x!=0, r)
        #assign equilibria to these positions
        sol_mat[k,inds] = real_sols[k]
    end
    #create column vector with subcommunity number
    subcomm_ind = repeat([j], n_sols)
    #create row vector and append to final data frame
    append_rows = hcat(subcomm_ind, sol_mat)
    global all_eq_mat = hcat(all_eq_mat, append_rows')
    print(all_eq_mat)
end #for (i)
#save them (by rows)
writedlm("../data/roots.csv", all_eq_mat')