# collision_integral.jl
# author: Davis Thuillier
# created: 10 July 2023

# include(joinpath("..", "src", "EEScattering2D.jl"))
using EEScattering2D

import StaticArrays: SVector
import LinearAlgebra: norm
import Dates
using Interpolations
using DelimitedFiles

function input_handling()
    if length(ARGS) != 5
        println("Insufficient number of arguments.")
        println("Format: collision_integral.jl [T::F64] [start::Int] [end::Int] [n_s::Int] [n_t::Int]")
        exit()
    end
end

function create_files(temperature::Float64, start_index::Int, end_index::Int, num_bins::Int, perp_num::Int)
    umklapp ? mit_umklapp = "umklapp" : mit_umklapp = "ohne_umklapp" 
    data_dir_stem = joinpath(@__DIR__, "..", "data", "$(band)_fermi_profile", "$(round(temperature, digits = 8))", "$(num_bins)_$(perp_num)")
    data_dir = data_dir_stem

    filenames = map( x -> "Γ" * x * "_$(row_dim)_$(round(temperature, digits = 8))_$(start_index)_to_$(end_index).csv", ["n","u"])
    
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
        # println("Data will be written to ", filenames[i])
    end
    
    return data_dir, filenames
end

function main() 
    data_dir, filenames = create_files(temperature, start_index, end_index, num_bins, perp_num)

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim)
    fv = Vector{SVector{2, Float64}}(undef, length(fs))
    FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  
    ds = FermiSurfaceMesh.get_ds(fs)

    ef::Float64 = 0.55 # Fermi energy in eV
    e0::Float64 = 55.26349406e6 # Vacuum permittivity in e^2 / eV / m
    c::Float64  = 12.68e-10 # z-axis dimension of unit cell in meters

    # Screened Parameters
    alpha = 1 / (ef * e0 * c) # Characteristic non-dimensionalized energy scale for interaction matrix element
    # prefactor = alpha^2 
    q_squared = alpha * get_dos(fs, fv) / (2pi)^2 # Thomas-Fermi screening wavevector squared 

    # Constant matrix element
    prefactor = (2pi)^2 / get_dos(fs, fv)

    open(joinpath(data_dir, "fermi_surface_$(length(fs)).csv"), "w") do file
        println(file, "kx,ky,vx,vy")
        writedlm(file, hcat(first.(fs), last.(fs), first.(fv), last.(fv)), ",")
    end

    interpolation_arclengths = FermiSurfaceMesh.get_arclengths(fs)
    perimeter = last(interpolation_arclengths)

    for i in start_index:end_index
        momenta, _ , _ , arclengths, _ = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, i, hamiltonian, temperature, prec)

        gamma = Matrix{Float64}(undef, size(momenta)[2] + 2, 3)
        
        FermiSurfaceIntegration.angular_distribution!(gamma, arclengths, perimeter, fs, momenta, hamiltonian, temperature, q_squared, umklapp = umklapp)

        # All values need to be multiplied by 2pi/(hbar / e_F)
        gamma = sortslices(gamma, dims = 1)
        gamma[begin, :] = gamma[end - 1, :] - [perimeter, 0.0, 0.0]
        gamma[end, :] = gamma[begin + 1, :] + [perimeter, 0.0, 0.0]

        itp_n = interpolate(gamma[:,1], gamma[:,2], FritschCarlsonMonotonicInterpolation())
        itp_u = interpolate(gamma[:,1], gamma[:,3], FritschCarlsonMonotonicInterpolation())

        # Write matrix data to file with momentum integrated out
        open(filenames[1], "a") do file
            writedlm(file, prefactor * sqrt(norm(fv[i]) * ds[i]) * transpose(itp_n.(interpolation_arclengths[1:end-1]) .* sqrt.(norm.(fv) .* ds)), ",")
        end

        open(filenames[2], "a") do file
            writedlm(file, prefactor * sqrt(norm(fv[i]) * ds[i]) * transpose(itp_u.(interpolation_arclengths[1:end-1]) .* sqrt.(norm.(fv) .* ds)), ",")
        end
        
    end
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))

input_handling()

const temperature = parse(Float64, ARGS[1])
const start_index = parse(Int,     ARGS[2])
const end_index   = parse(Int,     ARGS[3])
const num_bins    = parse(Int,     ARGS[4])
const perp_num    = parse(Int,     ARGS[5])

main()
