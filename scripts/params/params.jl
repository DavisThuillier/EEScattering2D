#### Parameters 
# const temperature::Float64 = 0.0006323
const row_dim::Int = 1200 # Number of points between 0 and 2pi excluding 2pi 

##### Constants #####
const prec::Float64   = 0.001 # Minimum value of F = f (1 - f) where f is the FD distribution to determine bounds of discretization
const tolerance::Float64  = 0.00001 
const max_iterations::Int = 10000
# const num_bins::Int       = 120 # Total number of angles to sample
# const perp_num::Int64     = 13 # The number of points along v_f for each angle; MUST BE AN EVEN NUMBER
const band::String        = "gamma"
const umklapp::Bool       = true 