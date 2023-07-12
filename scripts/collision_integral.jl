# collision_integral.jl
# author: Davis Thuillier
# created: 10 July 2023

using EEScattering2D

import StaticArrays: SVector
import LinearAlgebra: norm
import Dates
using Plots
using Interpolations
using DelimitedFiles
using ProgressBars

function create_files(start_index::Int, end_index::Int)
    data_dir_stem = joinpath(@__DIR__, "..", "data", band, "$(num_bins)_$(perp_num)", Dates.format(Dates.now(), "yyyy.mm.dd"))
    data_dir = data_dir_stem

    filenames = map( x -> "Γ" * x * "_$(row_dim)_$(round(temperature, digits = 4))_$(start_index)_to_$(end_index).csv", ["n","u", ""])
    
    # Check if any of the output files exist in the desired directory
    j = 0
    while sum(isfile.(joinpath.(data_dir, filenames))) > 0
        j += 1
        data_dir = data_dir_stem * "($(j))"
        !isdir(data_dir) && break
    end
    !isdir(data_dir) && mkpath(data_dir)

    for i in eachindex(filenames)
        filenames[i] = joinpath(data_dir, filenames[i])
        open(filenames[i], "w") do f
            # Clears file for writing
        end
        println("Data will be written to ", filenames[i])
    end
    
    return data_dir, filenames
end

function main(start_index::Int, end_index::Int) 
    data_dir, filenames = create_files(start_index, end_index)

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim)
    fv = Vector{SVector{2, Float64}}(undef, length(fs))
        FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  

    open(joinpath(data_dir, "fermi_surface_$(length(fs)).csv"), "w") do file
        println(file, "kx,ky,dh/dx,dh/dy")
        writedlm(file, hcat(first.(fs), last.(fs), first.(fv), last.(fv)), ",")
    end

    interpolation_arclengths = FermiSurfaceMesh.get_arclengths(fs)
    perimeter = FermiSurfaceMesh.get_perimeter(fs)

    for i in ProgressBar(start_index:end_index)
        @time momenta, dVs, variance, arclengths = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, i, hamiltonian, temperature, prec)

        plt = plot(first.(momenta), last.(momenta), seriestype= :scatter, markersize= 0.05, legend = false, markershape = :diamond, aspect_ratio = :equal, title = "T = $(round(temperature, digits = 4))", xlims = (-2.0,2.0), ylims = (-2.0,2.0))
        plot!(plt, first.(fs), last.(fs), color = :blue)
        display(plt)
        
        gamma = Matrix{Float64}(undef, length(arclengths) + 2, 3)
        
        FermiSurfaceIntegration.contracted_integral!(gamma, arclengths, perimeter, momenta, dVs, hamiltonian, variance, temperature)
        
        gamma = sortslices(gamma, dims = 1)
        gamma[begin, :] = gamma[end - 1, :] - [perimeter, 0.0, 0.0]
        gamma[end, :] = gamma[begin + 1, :] + [perimeter, 0.0, 0.0]

        itp_n = interpolate(gamma[:,1], gamma[:,2], FritschCarlsonMonotonicInterpolation())
        itp_u = interpolate(gamma[:,1], gamma[:,3], FritschCarlsonMonotonicInterpolation())

        # Write matrix data to file with momentum integrated out
        open(filenames[1], "a") do file
            writedlm(file, transpose(itp_n.(interpolation_arclengths)), ",")
        end

        open(filenames[2], "a") do file
            writedlm(file, transpose(itp_u.(interpolation_arclengths)), ",")
        end

        open(filenames[3], "a") do file
            writedlm(file, transpose(itp_n.(interpolation_arclengths) + itp_u.(interpolation_arclengths)), ",")
        end
        
    end
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))

main(1, div(row_dim, 8) + 1)