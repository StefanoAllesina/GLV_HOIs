using HomotopyContinuation #to solve system of polynomials
using LinearAlgebra #to compute dot products
using DelimitedFiles #to load and save files

function poly_i(x, r, A, B, i)
    #given parameters, construct glv_hoi system of polynomial equations
    r[i] + dot(A[i,:], x) + dot(x, B[:, :, i]*x)
end

function get_index(index, unit_history)
    #given the community number, and the richness previous communities
    #calculate the start index of such community in the parameter sets.
    from = Int16(sum(unit_history[1:index-1]) + 1)
    to = Int16(from + unit_history[index] - 1)
    [from, to]
end

#change working directory
cd("/Users/pablolechon/Desktop/phd/GLV_HOIs/code") #change this to your own path

#read data
#load growth rates, pairwise, three-way interactions, diversities of each comm
rs = readdlm("../data/rs.csv")
As = readdlm("../data/As.csv")
Bs = readdlm("../data/Bs.csv")
diversities = readdlm("../data/diversities.csv")
#total number of communities and biggest community
n_comms = size(diversities,1)
n_max = Int16(maximum(diversities))
#vector with column names of equilibria matrix
global all_eq_mat = Array{Float64}(undef, n_max)
global subcomm_ind = Vector{Float64}()
global n_spp_vec = Vector{Float64}()
#iterate through parameter sets (either communities or subcommunities)
for i in 1:n_comms
    #number of species in the ith (sub)community
    n_spp_i = Int16(diversities[i])
    #get indices for subsetting A and B
    ind_A = get_index(i, diversities)
    ind_B = get_index(i, diversities.^2)
    #get parameter set of ith (sub)community
    r = rs[i, 1:n_spp_i]
    A = As[ind_A[1]:ind_A[2], 1:n_spp_i]
    B = Bs[ind_B[1]:ind_B[2], 1:n_spp_i]
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
    sol_mat = zeros(n_sols, n_max)
    #loop through solutions
    for k in 1:n_sols
        #use growth rate vector to find present species 
        inds = findall(x->x!=0, r)
        #assign equilibria to these positions
        sol_mat[k,inds] = real_sols[k]
    end
    #create column vectors with subcommunity number and diversity of ith community
    global subcomm_ind = [subcomm_ind; repeat([i], n_sols)]
    global n_spp_vec = [n_spp_vec; repeat([n_spp_i], n_sols)]    
    global all_eq_mat = hcat(all_eq_mat, sol_mat')
end #for (i)
#delete first column (of zeros) of matrix 
all_eq_mat = all_eq_mat[:,2:end]
#add columns subcommunity index and community diversity on the left
all_eq_mat = hcat(subcomm_ind, n_spp_vec, all_eq_mat')
#save them (by rows)
writedlm("../data/roots.csv", all_eq_mat)
