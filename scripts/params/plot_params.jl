#### Plotting Parameters 
const temperatures = 0.01 * (2.0 .^ (0:0.5:2.5))
const num_angles::Int64 = 200
const p_num::Int64          = 20 # The number of points radial for each angle
const central_angle::Float64 = pi/4
const band::String          = "gamma"
const modes                  = 2:6
const mode_colors = [:blue :orange :green :pink :purple]

