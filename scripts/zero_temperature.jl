## author: Davis Thuillier
## about: Compute the momentum conserving pairs on the Fermi surface 

# Utility imports
import LinearAlgebra: norm, dot 
import Statistics: median, mean
import SpecialFunctions: erfinv, erf
using DelimitedFiles
using StaticArrays

using Base.Threads
using Interpolations

# Aesthetic imports
using ProgressBars 
using Plots

function hamiltonian(k::SVector{2, Float64})
    return - 2.0 * tx * cos(k[1]*pi) - 2.0 * ty * cos(k[2]*pi) - 4 * tp * cos(k[1]*pi) * cos(k[2]*pi)  - mu
end

fd(k::SVector{2, Float64}, T::Float64) = 1/(exp(hamiltonian(k)/T) + 1)

########################
#### Discretization ####
########################

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

function grad_hamiltonian(k::SVector{2,Float64})
    dp::Float64 = sqrt(eps(Float64))
    dh_x = hamiltonian(k + SVector{2}([dp,0])) - hamiltonian(k + SVector{2}([-dp, 0]))
    dh_y = hamiltonian(k + SVector{2}([0,dp])) - hamiltonian(k + SVector{2}([0, -dp]))
    return SVector{2}([dh_x, dh_y] / (2 * dp))
end

function fill_fermi_surface!(angles::Vector{Float64}, fermi_surface::Vector{SVector{2, Float64}})
    startpoint::SVector{2, Float64} = [0.0, 0.0]

    ## Compute roots of the Hamiltonian using the bisection method ##
    for i in eachindex(angles)
        n = SVector{2}([cos(angles[i]), sin(angles[i])])
        if - pi / 4 < angles[i] < pi / 4 || 3 * pi / 4 < angles[i] < 5 * pi / 4
            endpoint   = sqrt(1 + sin(angles[i])^2) * n
        else
            endpoint   = sqrt(1 + cos(angles[i])^2) * n
        end
        fermi_surface[i] = get_energy_root(startpoint, endpoint, 0.0)
    end

    return nothing
end

function fill_fermi_velocity!(fermi_surface::Vector{SVector{2, Float64}}, fermi_velocity::Vector{SVector{2, Float64}})
    for i in eachindex(fermi_surface)
        fermi_velocity[i] = grad_hamiltonian(fermi_surface[i])
    end
end

function count_pairs!(fermi_surface::Vector{SVector{2, Float64}}, num_pairs::Vector{Float64}, p1::SVector{2, Float64})
    P = SVector{2,Float64}
    K = SVector{2, Float64}
    norm_p1::Float64 = norm(p1)
    threshold = 0.01

    for i in eachindex(fermi_surface)
        P = p1 + fermi_surface[i] # fermi_surface[i] = p2
        K = p1 - fermi_surface[i] # fermi_surface[i]  = p1_prime
        norm_P = norm(P)
        norm_K = norm(K)
        for p1_prime in fermi_surface
            for p2_prime in fermi_surface
                (norm(P - p1_prime - p2_prime) / norm_p1 < threshold) && (num_pairs[i] -= norm(P - p1_prime - p2_prime) / norm_p1)
                (norm(K + p1_prime - p2_prime) / norm_p1 < threshold) &&(num_pairs[i] += 2 * norm(K + p1_prime - p2_prime) / norm_p1)
            end
        end
    end

    return nothing
end

function count_integral(fermi_surface::Vector{SVector{2, Float64}}, angles::Vector{Float64}, center_index::Int)
    pairs = zeros(Float64, length(angles))
    count_pairs!(fermi_surface, pairs, fermi_surface[center_index])
    pairs /= maximum(pairs)
    devs = Vector{SVector{2, Float64}}(undef, length(angles))
    for i in eachindex(devs)
        devs[i] = pairs[i] * SVector{2}([cos(angles[i]), sin(angles[i])])
    end

    plt = plot(angles, norm.(fermi_surface))
    plot!(plt, angles, norm.(fermi_surface .+ devs), proj = :polar)
    display(plt)
end

vector_atan(v::SVector{2,Float64}) = mod2pi(atan(v[2], v[1]))

rot(v::SVector{2, Float64}, theta::Float64) = SVector{2}([cos(theta) * v[1] - sin(theta) * v[2], sin(theta) * v[1] + cos(theta) * v[2]])

# function compute_weight(v1::SVector{2, Float64}, v2::SVector{2, Float64})
#     phi = mod2pi(atan( - (v1 - v2)[1], (v1 - v2)[2])) # Angle from 0 along which deviations in momentum conserve energy
#     psi = mod(phi - vector_atan(v1), pi)
#     n = rot(v1 / norm(v1), psi)

#     h_left = minimum([width / (2 * abs(cos(psi))), dp_plus / (2 * abs(sin(psi)))])
#     h_right = minimum([width / (2 * abs(cos(psi))), dp_minus / (2 * abs(sin(psi)))])
#     h = h_left + h_right

#     start = h_left * n
#     path = - h * n
#     temp_integral = 0.0
#     for k in 0:(p_num-1)
#         delta_p = start + (k / p_num) * path
#         temp_integral += (1 - fd(p1_prime + delta_p, temperature)) * (1 - fd(p2_prime - delta_p, temperature)) * (h / (p_num - 1))
#     end
# end

function main()
    include("params.jl")
    include(joinpath(@__DIR__,"bands", "tbm.jl"))

    N::Int = 200 # Index of pi/4
    total_num = 8*(N-1)
    dtheta = (pi/4) / (N-1)
    thetas = collect(range(0.0, 2 * pi, step = dtheta))
    pop!(thetas)

    fs = Vector{SVector{2, Float64}}(undef, length(thetas))
    fv = Vector{SVector{2, Float64}}(undef, length(thetas))
    fill_fermi_surface!(thetas, fs)
    fill_fermi_velocity!(fs, fv)

    # plot(first.(fs), last.(fs), aspect_ratio = 1.0)
    # plot!(first(fs[min_index]) .+ first.(fs), last(fs[min_index]) .+ last.(fs))
    integral = zeros(Float64, length(thetas))
    temp_integral::Float64 = 0.0

    p_num = 100
    mod_shift = SVector{2}([1.0,1.0])
        
    for i in 1:9:N
        p1 = fs[i] # Injection electron
        width = 2 * temperature
        for j in eachindex(fs)
            integral[j] = 0.0
            P = p1 + fs[j]
            K = p1 - fs[j]
            for l in eachindex(fs)
                p2_prime = mod.(P - fs[l] + mod_shift, 2.0) - mod_shift
                index2 = mod(Int(round(vector_atan( p2_prime ) / dtheta)), total_num) + 1 # Index of p2_prime
                process2 = norm(mod.(P - fs[l] - fs[index2] + mod_shift, 2.0) - mod_shift) / norm(p1) < sqrt(eps(Float64))

                p2_prime = mod.(K + fs[l] + mod_shift, 2.0) - mod_shift
                index34 = mod(Int(round(vector_atan(p2_prime) / dtheta)), total_num) + 1 # Index of p2_prime
                process34 = norm(mod.(K + fs[l] - fs[index34] + mod_shift, 2.0) - mod_shift) / norm(p1) < sqrt(eps(Float64))

                dp_plus = norm(fs[mod(l, total_num) + 1] - fs[l])
                dp_minus = norm(fs[l] - fs[mod(l - 2, total_num) + 1]) 

                if process2 
                    phi = mod2pi(atan( - (fv[l] - fv[index2])[1], (fv[l] - fv[index2])[2])) # Angle from 0 along which deviations in momentum conserve energy
                    psi = mod(phi - vector_atan(fv[l]), pi)
                    n = rot(fv[l] / norm(fv[l]), psi)

                    ### Distances from FS momentum to edge of discretization cell along energy conserving angle ### 
                    h_left = minimum([width / (2 * abs(cos(psi))), dp_plus / (2 * abs(sin(psi)))])
                    h_right = minimum([width / (2 * abs(cos(psi))), dp_minus / (2 * abs(sin(psi)))])
                    h = h_left + h_right
                
                    start = h_left * n
                    path = - h * n
                    temp_integral = 0.0
                    for k in 0:(p_num-1)
                        delta_p = start + (k / p_num) * path
                        temp_integral += (1 - fd(fs[l] + delta_p, temperature)) * (1 - fd(fs[index2] - delta_p, temperature)) * (h / (p_num - 1)) 
                    end
                    integral[j] += abs(cos(psi)) * (dp_minus + dp_plus) * fd(p1, temperature) * fd(fs[j], temperature) * temp_integral / temperature
                end

                if process34
                    phi = mod2pi(atan( - (fv[l] - fv[index34])[1], (fv[l] - fv[index34])[2])) # Angle from 0 along which deviations in momentum conserve energy
                    psi = mod(phi - vector_atan(fv[l]), pi)
                    n = rot(fv[l] / norm(fv[l]), psi)

                    ### Distances from FS momentum to edge of discretization cell along energy conserving angle ### 
                    h_left = minimum([width / (2 * abs(cos(psi))), dp_plus / (2 * abs(sin(psi)))])
                    h_right = minimum([width / (2 * abs(cos(psi))), dp_minus / (2 * abs(sin(psi)))])
                    h = h_left + h_right
                
                    start = h_left * n
                    path = - h * n
                    temp_integral = 0.0
                    for k in 0:(p_num-1)
                        delta_p = start + (k / p_num) * path
                        temp_integral += fd(fs[l] + delta_p, temperature) * (1 - fd(fs[index34] + delta_p, temperature)) * (h / (p_num - 1)) 
                    end
                    integral[j] += 2 * abs(cos(psi)) * (dp_minus + dp_plus) * fd(p1, temperature) * (1 - fd(fs[j], temperature)) * temp_integral / temperature
                end
            end
        end
        scale = 100 #0.02
        plt = plot(plot(thetas, norm.(fs), color = :black), proj = :polar)
        plot!(thetas, norm.(fs) .+ scale * integral, proj = :polar, ylims = (0,1.5)) 
        
        display( plt )

    end
    
    # open(joinpath(data_dir, "fermi_surface_$(length(fs)).csv"), "w") do file
    #     println(file, "kx,ky,dh/dx,dh/dy")
    #     writedlm(file, hcat(first.(fs), last.(fs), first.(fv), last.(fv)), ",")
    # end

end

@time main()