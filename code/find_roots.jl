using HomotopyContinuation #to solve system of polynomials
using LinearAlgebra #to compute dot products
using DelimitedFiles #to load and save files

function read_data(file_suffix)
    #read data given the number of parameter types and the suffix of files 
    rs = readdlm("../data/pars_1_"*file_suffix*".csv") #growth rates
    As = readdlm("../data/pars_2_"*file_suffix*".csv") #pairwise
    Bs = readdlm("../data/pars_3_"*file_suffix*".csv") #three-way 
    # diversities of each comm
    diversities = readdlm("../data/pars_4_"*file_suffix*".csv") 
    #species pool vector
    pools = readdlm("../data/pars_5_"*file_suffix*".csv") 
    rs, As, Bs, diversities, pools
end

function get_index(index, unit_history)
    #given the community number, and the richness previous communities
    #calculate the start index of such community in the parameter sets.
    from = Int16(sum(unit_history[1:index-1]) + 1)
    to = Int16(from + unit_history[index] - 1)
    [from, to]
end

function poly_i(x, r, A, B, i)
    #given parameters, construct glv_hoi system of polynomial equations
    r[i] + dot(A[i,:], x) + dot(x, B[:, :, i]*x)
end

function solve_system(equations, presence)
    #instantiate polynomial system 
    F = System(equations)
    #initialize variables
    n_sols = Array{Int64}(undef)
    sol_mat = Array{Float64}(undef, 0)
    try
        #try to solve polynomial system
        result = solve(F) #might throw error
        #get number of real solutions
        real_sols = real_solutions(result)
        n_sols = size(real_sols, 1)
        #preallocate n_sols rows of zeros to store solutions
        sol_mat = zeros(n_sols, n_max)
        #loop through solutions
        for k in 1:n_sols
            #use growth rate vector to find present species 
            inds = findall(x->x!=0, presence)
            #assign solutions to these positions
            sol_mat[k,inds] = real_sols[k]
        end
    catch
       #if error when solving the system, flag it as a row of -1s
       n_sols = 1
       sol_mat = -1*ones(Int, (n_sols, n_max))
       println("couldn't find solutions of system")
    end
    n_sols, sol_mat
end

#parse of output parameter file name suffix spedifying parameter constraints
output_name = ARGS[1]
#read data
rs, As, Bs, diversities, pools = read_data(output_name)
#total number of communities and biggest community
n_comms = size(diversities,1)
n_max = Int16(maximum(diversities))
#vector with column names of equilibria matrix
global all_eq_mat = Array{Float64}(undef, n_max)
#vector with community diversity
global n_spp_vec = Vector{Float64}()
#vector with solution number
global sol_id = Vector{Int64}()
#vector labeling species pool (common to a community and all subcommunities)
global pool_id_vec = Vector{Int64}()
#iterate through parameter sets of each (sub)community
for i in 1:n_comms
    print("Community: ")
    print(i)
    print(", Diversity: ")
    println(diversities[i])
    #number of species in the ith (sub)community
    n_spp_i = Int16(diversities[i])
    #assign pool id for this iteration
    pool_id_i = Int16(pools[i])
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
    #solve system
    n_sols, sol_mat = solve_system(equations, r)
    #create column vectors with solution id, diversity of ith community, 
    #and species pool id	
    global sol_id = [sol_id; repeat([i], n_sols)]
    global n_spp_vec = [n_spp_vec; repeat([n_spp_i], n_sols)] 
    global pool_id_vec = [pool_id_vec; repeat([pool_id_i], n_sols)] 
    #store solutions in the global matrix   
    global all_eq_mat = hcat(all_eq_mat, sol_mat')
end #for (i)
#delete first column (of zeros) of matrix 
all_eq_mat = all_eq_mat[:,2:end]
#prepend columns iteration counter and community diversity
all_eq_mat = hcat(sol_id, n_spp_vec, pool_id_vec, all_eq_mat')
#save them by rows
writedlm("../data/roots_"*output_name*".csv", all_eq_mat)
