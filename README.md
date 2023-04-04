# GLV_HOIs
Generalized Lotka-Volterra model with Higher-Order Interactions
the following files numerically compute all interior and boundary equilibria of a GLV system with HOIs

1. get_all_equilibria.r  is the template code
2. build_par_subcomms.r  generates the 2^n -1 parameter sets of all sub-communities of the given species poll
3. find_roots.jl  finds all equilibria of all sub-systems using Homotopy Continuation in one go.
