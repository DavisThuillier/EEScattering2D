## author: Davis Thuillier
## about: Compute the Boltzmann collision integral for ee scattering in the basis of delta functions of angle from the injected particle momentum. Output is a CSV file for each temperature specified in params.jl

# Utility imports
import LinearAlgebra: norm, dot 
import Statistics: median, mean
import SpecialFunctions: erfinv, erf
using DelimitedFiles
using StaticArrays

using Interpolations

# Aesthetic imports
using ProgressBars 
using Plots
# using LaTeXStrings

function hamiltonian(k::SVector{2, Float64})
    #return 0.5 * (k[1]^2 + k[2]^2) + mu # Free electron gas
    return - 2.0 * tx * cos(k[1]*pi) - 2.0 * ty * cos(k[2]*pi) - 4 * tp * cos(k[1]*pi) * cos(k[2]*pi)  - mu
end

fd(k::SVector{2, Float64}, T::Float64) = 1/(exp(hamiltonian(k)/T) + 1)

########################
#### Discretization ####
########################

collinear_width(T::Float64) = max_width * (1 - (1 - min_width/max_width) * exp( - (T - 0.0025) * width_slope))

function gaussian_density(x::Float64, sigma::Float64, amplitude::Float64, limit_ratio::Float64)
    res = amplitude * exp( - (x / sigma)^2)
    if res <= limit_ratio * amplitude
        return limit_ratio * amplitude
    else 
        return res
    end 
end 

"Amplitude for density function normalized for N points between 0 point and endpoint."
function get_amp(N::Float64, sigma::Float64, ratio::Float64, endpoint::Float64)
    return N / (sqrt(pi) * sigma * erf(log(1/ratio)) / 2 + (endpoint - sigma * log(1/ratio)) * ratio)
end

function get_energy_root(startpoint::SVector{2, Float64}, endpoint::SVector{2, Float64}, energy::Float64 = 0.0)
    delta_E::Float64 = 0.0
    n = endpoint - startpoint
    n = n / norm(n)

    j::Int = 0
    while j < max_iterations
        midpoint = (startpoint + endpoint) / 2
        delta_E = hamiltonian(midpoint) - energy
        if norm(endpoint - startpoint) / 2 < tolerance
            return midpoint
        end
        
        if sign(delta_E) == sign(hamiltonian(startpoint))
            startpoint = midpoint
        else
            endpoint   = midpoint
        end

        j += 1
    end

    return (startpoint + endpoint) / 2
end

function fill_fermi_surface!(angles::Vector{Float64}, fermi_surface::Vector{SVector{2, Float64}})
    startpoint::SVector{2, Float64} = [0.0, 0.0]

    ## Compute roots of the Hamiltonian using the bisection method ##
    for i in eachindex(angles)
        n = SVector{2}([cos(angles[i]), sin(angles[i])])
        endpoint = 2.0 * n 
        if - pi / 4 < angles[i] < pi / 4 || 3 * pi / 4 < angles[i] < 5 * pi / 4
            endpoint   = sqrt(1 + sin(angles[i])^2) * n
        else
            endpoint   = sqrt(1 + cos(angles[i])^2) * n
        end
        fermi_surface[i] = get_energy_root(startpoint, endpoint, 0.0)
    end

    return nothing
end

function grad_hamiltonian(k::SVector{2,Float64})
    dp::Float64 = sqrt(eps(Float64))
    dh_x = hamiltonian(k + SVector{2}([dp,0])) - hamiltonian(k + SVector{2}([-dp, 0]))
    dh_y = hamiltonian(k + SVector{2}([0,dp])) - hamiltonian(k + SVector{2}([0, -dp]))
    return SVector{2}([dh_x, dh_y] / (2 * dp))
end

function fill_fermi_velocity!(fermi_surface::Vector{SVector{2, Float64}}, fermi_velocity::Vector{SVector{2, Float64}})
    for i in eachindex(fermi_surface)
        fermi_velocity[i] = grad_hamiltonian(fermi_surface[i])
    end
    return nothing
end

function get_arclengths(curve::Vector{SVector{2, Float64}})
    arc_lengths = Vector{Float64}(undef, length(curve))
    arc_lengths[1] = 0.0

    for i in 2:length(curve)
        arc_lengths[i] = arc_lengths[i - 1] + norm(curve[i] - curve[i - 1])
    end

    return arc_lengths
end

function get_angles(fermi_surface::Vector{SVector{2, Float64}}, injection_angle_index::Int, num_angles::Int, sigma::Float64)
    centered_fs = circshift(fermi_surface, 1 - injection_angle_index) # Make the injection angle the first element of the array
    centered_fs_r = circshift(reverse(centered_fs), 1)
    t = get_arclengths(centered_fs) # Compute arc_lengths from the central angle
    t_r = get_arclengths(centered_fs_r)
    sigma = sigma * t[end] / (2*pi) # Scale width of the peak by total arclength compared to unit circle
    central_angle = mod(atan(centered_fs[1][2], centered_fs[1][1]), 2 *pi)
    angles::Vector{Float64} = [central_angle]
    angles_r::Vector{Float64} = []
    limit_ratio = 0.2

    amp = get_amp(num_angles / 4, sigma, limit_ratio, t[end] / 4) 
    
    s::Float64     = 0.0
    i::Int         = 1
    theta::Float64 = central_angle
    k::Vector{Float64} = [0.0,0.0]

    extend_domain = (central_angle + pi >= 2*pi)
    while (theta < pi/2 + central_angle) && i < length(t)
        s += 1 / gaussian_density(s, sigma, amp, limit_ratio)
        while s > t[i]
            i += 1
        end
        k = centered_fs[i - 1] + ( (s - t[i - 1]) / (t[i] - t[i-1])) * (centered_fs[i] - centered_fs[i - 1])
        theta = mod(atan(k[2], k[1]), 2*pi)
        extend_domain && ((theta <= pi) && (theta += 2*pi)) 
        (theta < pi/2 + central_angle) && push!(angles, theta)
    end

    s = 0.0
    i = 1
    theta = central_angle
    extend_domain = (central_angle - pi <= 0)
    while theta > (central_angle - pi/2) && i < length(t_r)
        s += 1 / gaussian_density(s, sigma, amp, limit_ratio)
        while s > t_r[i]
            i += 1
        end
        k = centered_fs_r[i - 1] + ( (s - t_r[i - 1]) / (t_r[i] - t_r[i-1])) * (centered_fs_r[i] - centered_fs_r[i - 1])
        theta = mod(atan(k[2], k[1]), 2*pi)
        extend_domain && ((theta > pi) && (theta -= 2*pi))

        theta > (central_angle - pi/2) && push!(angles_r, theta)
    end

    (length(angles_r) > num_angles / 2) && err("Discretization error: too many angles.")
    (length(angles) > num_angles / 2) && err("Discretization error: too many angles.")

    return vcat(angles, pi .+ reverse(angles_r), pi .+ angles, reverse(angles_r))
end

function get_k_bound(e_bound::Float64, fs_k::SVector{2, Float64}, velocity::SVector{2,Float64})
    n = velocity / norm(velocity)
    step = (e_bound / norm(velocity)) / 2

    i::Int = 0
    i_limit::Int = abs(div(1, step)) # Number of step corresponding to half-width of Brillouin zone

    endpoint = fs_k
    startpoint = fs_k
    
    if e_bound > 0
        while i < i_limit && hamiltonian(endpoint) < e_bound
            endpoint += step * n
            i += 1
        end
    else
        while i < i_limit && hamiltonian(startpoint) > e_bound
            startpoint += step * n
            i += 1
        end
    end

    j::Int = 0
    while j < max_iterations
        midpoint = (startpoint + endpoint) / 2
        delta_E = hamiltonian(midpoint) - e_bound
        norm(endpoint - startpoint) / 2 < tolerance && break
        
        if sign(delta_E) == sign(hamiltonian(startpoint) - e_bound)
            startpoint = midpoint
        else
            endpoint   = midpoint
        end
        j += 1
    end

    k_bound = (startpoint + endpoint) / 2

    #Check if k_bound lies outside the Brillouin zone
    if abs(k_bound[1]) > 1
        k_bound = fs_k + ( (sign(k_bound[1]) - fs_k[1]) / n[1]) * n
    elseif abs(k_bound[2]) > 1
        k_bound = fs_k + ( (sign(k_bound[2]) - fs_k[2]) / n[2]) * n
    end
    return k_bound
end

function temperature_broaden(fermi_surface::Vector{SVector{2, Float64}}, fermi_velocity::Vector{SVector{2, Float64}}, T::Float64)
    e_max::Float64 = 2  * T * acosh(1 / (2 * sqrt(prec)))

    momenta = Matrix{SVector{2, Float64}}(undef, p_num, length(fermi_surface))

    for i in eachindex(fermi_surface)
        # p_max = get_k_bound(e_max, n, fermi_surface[i], step)
        p_min = get_k_bound(-e_max, fermi_surface[i], fermi_velocity[i])
        p_max = get_k_bound(e_max, fermi_surface[i], fermi_velocity[i])

        for j in 1:p_num # 2 additional points correspond to boundary contours of ring
            @inbounds momenta[j, i] = p_min + (p_max - p_min) * (j - 1) / (p_num - 1)
        end
    end

    return momenta
end

area(v::SVector{2, Float64}, u::SVector{2, Float64}) = abs(v[1] * u[2] - v[2] * u[1]) # Area in the plane spanned by 2-vectors

# Stores the indices in the momentum array associated to a given angle
struct Bin
    angle::Float64 
    indices::Vector{Tuple{Int, Int}}
end

function discretize(fermi_surface::Vector{SVector{2, Float64}}, injection_angle_index::Int, T::Float64)
    colin_width::Float64 = collinear_width(T)
    angles = mod2pi.( get_angles(fermi_surface, injection_angle_index, num_angles,colin_width / (2 * sqrt(2 * log(2)))) ) # Get new angles on which to discretize

    new_fs = Vector{SVector{2, Float64}}(undef, length(angles))
    new_fv = Vector{SVector{2, Float64}}(undef, length(angles))
    fill_fermi_surface!(angles, new_fs)
    fill_fermi_velocity!(new_fs, new_fv)

    dV      = zeros(Float64, p_num, length(new_fs))
    momenta = temperature_broaden(new_fs, new_fv, T)
    
    angle_bins = Vector{Bin}(undef, length(angles))
    for i in eachindex(angles)
        angle_bins[i] = Bin(angles[i], [])
    end

    for i in 2:(size(momenta)[1] - 1) # Ignore boundary contours
        for j in 1:size(momenta)[2]
            dV[i,j] = area((momenta[i + 1, j] - momenta[i - 1, j])/ 2, (momenta[i, mod(j, length(new_fs)) + 1] - momenta[i, mod(j - 2, length(new_fs)) + 1]) / 2)

            theta = mod2pi(atan(momenta[i,j][2], momenta[i,j][1])) 
            min_index = argmin( abs.( - theta .+ angles))
            push!(angle_bins[min_index].indices, (i,j))
        end
    end

    for bin in angle_bins
        sort!(bin.indices, by = (x) -> norm(momenta[x[1], x[2]]))
    end

    variance = median(dV) / 4 # Approximate minimal width squared for delta function normalization

    return momenta, dV, angle_bins, variance
end

#####################
#### Integration ####
#####################

fd_normalization(k::SVector{2, Float64}, T::Float64) = 4 * cosh(hamiltonian(k) / (2 * T))

gaussian_delta(deviation::Float64, sigma_squared::Float64) = exp( - deviation^2 / (2 * sigma_squared)) / (sigma_squared * sqrt(2 * pi))

function collision_integral_riemann(p1::SVector{2, Float64}, k::SVector{2, Float64}, momenta::Matrix{SVector{2, Float64}}, dVs::Matrix{Float64}, sigma_squared::Float64, T::Float64)
    I2_n::Float64 = 0.0 
    I2_u::Float64 = 0.0
    I34_n::Float64 = 0.0
    I34_u::Float64 = 0.0

    P = p1 + k # Total momentum of scatterers
    Q = p1 - k # Momentum difference for p1 -> p1'

    E2::Float64 = hamiltonian(p1) + hamiltonian(k)
    E34::Float64 = hamiltonian(p1) - hamiltonian(k)

    mod_shift = SVector{2}([1.0,1.0]) # For shifting the wavevector over before taking the modulus with respect to the first Brillouin Zone

    for i in 2:(size(momenta)[1] - 1)
        for j in 1:size(momenta)[2]
            @inbounds p1_prime = momenta[i,j]
            p2_prime = mod.(mod_shift + P - p1_prime, 2.0) - mod_shift
            #p2_prime = P - p1_prime
            if ( abs((P - p1_prime)[1] > 1) || abs((P - p1_prime)[1] > 1) ) # Detect whether this is an umklapp event
                I2_u += fd(k, T) * (1 - fd(p1_prime, T)) * (1 - fd(p2_prime, T)) * dVs[i] * gaussian_delta(E2 - hamiltonian(p1_prime) - hamiltonian(p2_prime), sigma_squared) 
            else
                I2_n += fd(k, T) * (1 - fd(p1_prime, T)) * (1 - fd(p2_prime, T)) * dVs[i] * gaussian_delta(E2 - hamiltonian(p1_prime) - hamiltonian(p2_prime), sigma_squared) 
            end

            @inbounds p2 = momenta[i,j]
            p2_prime = mod.(mod_shift + Q + p2, 2) - mod_shift
            #p2_prime = Q + p2
            if ( abs((Q + p2)[1] > 1) || abs((Q + p2)[1] > 1) ) # Detect whether this is an umklapp event
                I34_n += 2 * fd(p2, T) * (1 - fd(k, T)) * (1 - fd(p2_prime, T)) * dVs[i] * gaussian_delta((E34 + hamiltonian(p2) - hamiltonian(p2_prime)), sigma_squared)
            else
                I34_u += 2 * fd(p2, T) * (1 - fd(k, T)) * (1 - fd(p2_prime, T)) * dVs[i] * gaussian_delta((E34 + hamiltonian(p2) - hamiltonian(p2_prime)), sigma_squared)
            end
        end
    end
    
    return [(- I2_n + I34_n) * fd(p1, T),  (- I2_u + I34_u) * fd(p1, T)] # fd(p1, T) is a constant in each integral and thus removed
end

function get_dp(momenta::Matrix{SVector{2,Float64}}, bin::Bin)
    # Returns the radial momentum differential assuming sorted momenta in the bin
    bin_momenta = map((x) -> momenta[x[1], x[2]], bin.indices)
    dp = Vector{Float64}(undef, sizeof(bin_momenta))

    for i in eachindex(bin_momenta)
        if i < length(bin_momenta)
            dp[i] = norm(bin_momenta[i + 1]) - norm(bin_momenta[i])
        else 
            dp[i] = dp[i - 1]           
        end
    end

    return dp
end

"Compute Boltzmann collision integral between each binned angle and the injection angle."
function angular_integral!(momenta::Matrix{SVector{2, Float64}}, dVs::Matrix{Float64}, bins::Vector{Bin}, variance::Float64, T::Float64, reduced_mat::Matrix{Float64})
    integral::SVector{2,Float64} = [0.0,0.0]

    central_momenta::Vector{SVector{2, Float64}} = map((x) -> momenta[x[1], x[2]], bins[1].indices)
    central_dp::Vector{Float64} = get_dp(momenta, bins[1])
    
    for i in ProgressBar(eachindex(bins))
        dp = get_dp(momenta, bins[i])
        integral = SVector{2}([0.0, 0.0])
        for m in eachindex(bins[i].indices)
            coord = bins[i].indices[m]
            k = momenta[coord[1], coord[2]]
            for j in eachindex(central_momenta)
                p1 = central_momenta[j]
                loss_terms = collision_integral_riemann(p1, k, momenta, dVs, variance, T) * fd_normalization(p1, T)

                integral += loss_terms * (norm(p1) * norm(k) * dp[m] * central_dp[j]) / T # Divide by T due to delta functions in the integration
            end   
        end
        reduced_mat[i, 1] = bins[i].angle
        reduced_mat[i, 2] = integral[1] # Normal scattering contribution
        reduced_mat[i, 3] = integral[2] # Umklapp scattering contribution
    end

    ### Padding rows ###
    reduced_mat[end - 1, :] = [-6*pi, 0.0, 0.0] # Will be sorted as the first element
    reduced_mat[end, :] = [6*pi, 0.0, 0.0] # Will be sorted as the final element
    return nothing
end

macro Name(variable)
    string(variable)
end

function single_run(parameter_file::String, central_angle::Float64; display_mesh::Bool = false)
    include(parameter_file)
    include(joinpath(@__DIR__,"bands", "$(band).jl"))

    dtheta = (pi/4) / (row_dim-1)
    thetas = collect(range(0.0, 2 * pi, step = dtheta))
    pop!(thetas)

    data_dir = joinpath(@__DIR__, "data", band) 
    data_dir = joinpath(data_dir, "$(num_angles)_$(p_num)")
    !isdir(data_dir) && mkpath(data_dir)

    fs = Vector{SVector{2, Float64}}(undef, length(thetas))
    fv = Vector{SVector{2, Float64}}(undef, length(thetas))
    fill_fermi_surface!(collect(thetas), fs)
    fill_fermi_velocity!(fs, fv)
    open(joinpath(data_dir, "fermi_surface_$(length(fs)).csv"), "w") do file
        println(file, "kx,ky,dh/dx,dh/dy")
        writedlm(file, hcat(first.(fs), last.(fs), first.(fv), last.(fv)), ",")
    end

    min_index = argmin(abs.(-central_angle .+ thetas))
    single_run_dir = joinpath(data_dir, "theta_$(round(thetas[min_index], digits=4))")
    !isdir(single_run_dir) && mkpath(single_run_dir)

    grid_vals = Base._linspace(-1.0, 1.0, 200)
    grid = collect(Iterators.product(grid_vals, grid_vals))
    
    for t in temperatures
        @time momenta, dVs, bins, variance = discretize(fs, min_index, t)

        if display_mesh
            plt = heatmap(grid_vals, grid_vals, hamiltonian.(SVector{2}.(grid)), colormap = :nuuk, colorbar_title = latexstring("\$\\epsilon(\\vec{k}) / t_x\$"))
            plot!(plt, first.(momenta), last.(momenta), seriestype= :scatter, markersize= 0.05, legend = false, markershape = :diamond, aspect_ratio = :equal, title = "T = $(round(t, digits = 4))", xlims = (-1,1), ylims = (-1,1))
            plot!(plt, first.(fs), last.(fs), color = :blue)
            display(plt)
        end

        I_angular = Matrix{Float64}(undef, length(bins), 3) # Columns: theta, I_n (normal scattering), I_u (umklapp scattering)
        
        @time angular_integral!(momenta, dVs, bins, variance, t, I_angular)

        open(joinpath(single_run_dir, "T_$(round(t,digits = 4)).csv"), "w") do file
            println(file, "theta,I_n,I_u")
            writedlm(file, I_angular, ",")
        end
    end 
end

function main(parameter_file::String, start_index::Int, end_index::Int)
    include(parameter_file)
    include(joinpath(@__DIR__,"bands", "$(band).jl"))

    dtheta = (pi/4) / (row_dim-1)
    thetas = collect(range(0.0, 2 * pi, step = dtheta))

    data_dir = joinpath(@__DIR__, "data", band) 
    data_dir = joinpath(data_dir, "$(num_angles)_$(p_num)")
    !isdir(data_dir) && mkpath(data_dir)

    fs = Vector{SVector{2, Float64}}(undef, length(thetas))
    fv = Vector{SVector{2, Float64}}(undef, length(thetas))
    fill_fermi_surface!(collect(thetas), fs)
    fill_fermi_velocity!(fs, fv)
    open(joinpath(data_dir, "fermi_surface_$(length(fs)).csv"), "w") do file
        println(file, "kx,ky,dh/dx,dh/dy")
        writedlm(file, hcat(first.(fs), last.(fs), first.(fv), last.(fv)), ",")
    end

    filenames = map( x -> joinpath(data_dir, "Γ" * x * "_$(length(thetas))_$(round(temperature, digits = 4))_$(start_index)_to_$(end_index).csv"), ["n","u", ""])
    for file in filenames
        open(file, "w") do f
        end
        println("Data will be written to ", file)
    end
    
    #(1 == 1) && return nothing
    I_interp_n = Vector{Float64}(undef, length(thetas))
    I_interp_n = Vector{Float64}(undef, length(thetas))
    for i in start_index:end_index
        momenta, dVs, bins, variance = discretize(fs, i, temperature)

        plt = plot(first.(momenta), last.(momenta), seriestype= :scatter, markersize= 0.05, legend = false, markershape = :diamond, aspect_ratio = :equal, title = "T = $(round(temperature, digits = 4))", xlims = (-2.0,2.0), ylims = (-2.0,2.0))
        plot!(plt, first.(fs), last.(fs), color = :blue)
        display(plt)
        
        I_angular = Matrix{Float64}(undef, length(bins) + 2, 3)
        
        angular_integral!(momenta, dVs, bins, variance, temperature, I_angular)

        I_angular = sortslices(I_angular, dims = 1)
        I_angular[begin, :] = I_angular[end - 1, :] - [2*pi, 0.0, 0.0]
        I_angular[end, :] = I_angular[begin + 1, :] + [2*pi, 0.0, 0.0]
        
        itp_n = interpolate(I_angular[:,1], I_angular[:,2], FritschCarlsonMonotonicInterpolation())
        itp_u = interpolate(I_angular[:,1], I_angular[:,3], FritschCarlsonMonotonicInterpolation())

        # Write matrix data to file with momentum integrated out
        open(filenames[1], "a") do file
            writedlm(file, transpose(itp_n.(thetas)), ",")
        end

        open(filenames[2], "a") do file
            writedlm(file, transpose(itp_u.(thetas)), ",")
        end

        open(filenames[3], "a") do file
            writedlm(file, transpose(itp_n.(thetas) + itp_u.(thetas)), ",")
        end

        #println("Data recorded to ", filename)
        
    end
end

# function error_handling()
#     if length(ARGS) != 3
#         println("Input structure: julia main.jl [params_file.jl] [start_index] [end_index]")
#         exit()
#     end

#     if !isfile(ARGS[1])
#         println("Parameter file does not exist.")
#         exit()
#     end

#     if isodd(p_num) || !isinteger(p_num) || p_num < 1 
#         error("$(@Name p_num) must be an even, positive integer.")
#     end
# end

#main(ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]))
main("params.jl", 1, 101)
#single_run("params.jl", pi/6; display_mesh = true)
#single_run("params.jl", pi/6; display_mesh = true)