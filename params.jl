#### Parameters 
#const temperatures::Vector{Float64} = 0.01 * (2.0 .^ (0:0.5:2.5))
const temperature::Float64 = 0.04
const row_dim::Int = 101 # Number of points between 0 and pi/4 excluding pi/4 

##### Constants #####
const prec::Float64   = 0.001 # Minimum value of F = f (1 - f) where f is the FD distribution to determine bounds of discretization
const tolerance::Float64   = 0.00001 
const max_iterations::Int = 10000
const num_angles::Int       = 200 # Total number of angles to sample
const p_num::Int64          = 30 # The number of points along v_f for each angle; MUST BE AN EVEN NUMBER
const band::String          = "gamma"

## Collinear width parameters ##
const max_width::Float64 = 0.4 # Maximum width of collinear region in radians
const min_width::Float64 = 0.2 
const width_slope::Float64 = 1.0 # Initial slope of exponential plateau function
#####################