# collision_integral.jl
# author: Davis Thuillier
# created: 10 July 2023

using EEScattering2D

import StaticArrays: SVector
using Plots
using Interpolations

function main(start_index::Int, end_index::Int)
    dtheta = (pi/4) / (row_dim-1)
    thetas = collect(range(0.0, 2 * pi, step = dtheta))

    data_dir = joinpath(@__DIR__, "..", "data", band, "$(num_bins)_$(perp_num)")
    !isdir(data_dir) && mkpath(data_dir)

    fs, _ = uniform_fermi_surface(thetas, hamiltonian, data_dir, write = true)

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
        @time momenta, dVs, bins, variance = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, i, hamiltonian, temperature, prec)

        plt = plot(first.(momenta), last.(momenta), seriestype= :scatter, markersize= 0.05, legend = false, markershape = :diamond, aspect_ratio = :equal, title = "T = $(round(temperature, digits = 4))", xlims = (-2.0,2.0), ylims = (-2.0,2.0))
        plot!(plt, first.(fs), last.(fs), color = :blue)
        display(plt)
        
        # I_angular = Matrix{Float64}(undef, length(bins) + 2, 3)
        
        # angular_integral!(momenta, dVs, bins, variance, temperature, I_angular)

        # I_angular = sortslices(I_angular, dims = 1)
        # I_angular[begin, :] = I_angular[end - 1, :] - [2*pi, 0.0, 0.0]
        # I_angular[end, :] = I_angular[begin + 1, :] + [2*pi, 0.0, 0.0]
        
        # itp_n = interpolate(I_angular[:,1], I_angular[:,2], FritschCarlsonMonotonicInterpolation())
        # itp_u = interpolate(I_angular[:,1], I_angular[:,3], FritschCarlsonMonotonicInterpolation())

        # # Write matrix data to file with momentum integrated out
        # open(filenames[1], "a") do file
        #     writedlm(file, transpose(itp_n.(thetas)), ",")
        # end

        # open(filenames[2], "a") do file
        #     writedlm(file, transpose(itp_u.(thetas)), ",")
        # end

        # open(filenames[3], "a") do file
        #     writedlm(file, transpose(itp_n.(thetas) + itp_u.(thetas)), ",")
        # end
        
    end
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))

main(1, 5)