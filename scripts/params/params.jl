#### Parameters 
const temperature::Float64 = 0.04
const row_dim::Int = 101 # Number of points between 0 and pi/4 excluding pi/4 

##### Constants #####
const prec::Float64   = 0.001 # Minimum value of F = f (1 - f) where f is the FD distribution to determine bounds of discretization
const tolerance::Float64   = 0.00001 
const max_iterations::Int = 10000
const num_bins::Int       = 50 # Total number of angles to sample
const perp_num::Int64          = 10 # The number of points along v_f for each angle; MUST BE AN EVEN NUMBER
const band::String          = "gamma"