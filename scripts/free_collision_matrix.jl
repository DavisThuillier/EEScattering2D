using EEScattering2D

import StaticArrays: SVector
import LinearAlgebra: norm, dot
using Interpolations
using DelimitedFiles
using ProgressBars

function input_handling()
    if length(ARGS) != 6
        println("Insufficient number of arguments.")
        println("Format: com_collision_matrix.jl [band::Str] [T::F64] [n_s::Int] [n_t::Int] [n_ξ1::Int] [dim::Int]")
        exit()
    end
end

function create_files(temperature::Float64, num_bins::Int, perp_num::Int, num_samples::Int, row_dim::Int)
    data_dir_stem = joinpath(@__DIR__, "..", "data", "$(band)_band", "$(round(temperature, digits = 8))", "$(num_bins)x$(perp_num)_$(num_samples)")
    data_dir = data_dir_stem

    filenames = map(x -> "Γ_$(x)_$(row_dim)_$(round(temperature, digits = 8)).csv", ["s"]) 

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
        open(filenames[i], "w") do file
            # Clear output file for writing 
        end
    end
    
    
    return data_dir, filenames
end

function main()
    data_dir, outfiles = create_files(temperature, num_bins, perp_num, sym_num, row_dim)

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim, bz = bz) # Uniform Fermi Surface
    fv = Vector{SVector{2, Float64}}(undef, length(fs))
    FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian) 
    ds = FermiSurfaceMesh.get_ds(fs) 
    arclengths = FermiSurfaceMesh.get_arclengths(fs)
    perimeter = last(arclengths)
    pop!(arclengths) # Remove the final element so that arclengths can be used as the interpolation mesh

    ############################
    ## Thomas-Fermi Screening ##
    ############################

    ef::Float64 = 0.55 # Fermi energy in eV
    e0::Float64 = 55.26349406e6 # Vacuum permittivity in e^2 / eV / m
    c::Float64  = 12.68e-10 # z-axis dimension of unit cell in meters

    # Screened Parameters
    alpha = 1 / (ef * e0 * c) # Characteristic non-dimensionalized energy scale for interaction matrix element
    q_squared = alpha * get_dos(fs, fv) / (2pi)^2 # Thomas-Fermi screening wavevector squared 
    prefactor = (2pi)^2 / get_dos(fs, fv) # For a constant interaction matrix element

    # Write referenced FS and corresponding Fermi velocity to file
    open(joinpath(data_dir, "fermi_surface_$(length(fs)).csv"), "w") do file
        println(file, "kx,ky,vx,vy")
        writedlm(file, hcat(first.(fs), last.(fs), first.(fv), last.(fv)), ",")
    end
    amp_ratio = 8.0

    s1 = 0.0
    for i in 
        Γ_s[i,j] = prefactor * sum(FermiSurfaceIntegration.contracted_integral(s1, s2, fs, num_bins, perp_num, hamiltonian, temperature, q_squared, umklapp = umklapp, bz = false))

    
    s1::Float64 = 0.0
    s2::Float64 = 0.0
    Γ_ξ = Matrix{Float64}(undef, asym_num, sym_num) 
    for (j, ξ1) in ProgressBar(enumerate(sym_mesh))
        for (i, ξ2) in enumerate(asym_mesh)
            s1 = mod((ξ1 + ξ2) / sqrt(2), perimeter)
            s2 = mod((ξ1 - ξ2) / sqrt(2), perimeter)
            
            Γ_ξ[i,j] = prefactor * sum(FermiSurfaceIntegration.contracted_integral(s1, s2, fs, num_bins, perp_num, hamiltonian, temperature, q_squared, umklapp = umklapp, bz = bz)) # Sums the normal and umklapp contributions
        end
    end

    open(outfiles[1], "w") do file
        writedlm(file, Γ_ξ, ",")
    end

    ###################
    ## Interpolation ##
    ###################
    itp = interpolate((asym_mesh, sym_mesh), Γ_ξ, Gridded(Linear()))
    etp = extrapolate(itp, Flat())

    # Collision matrix in (s1,s2) coordinates
    # s1 ↦ row 
    # s2 ↦ column
    Γ_s = Matrix{Float64}(undef, row_dim, row_dim)
    sym_factor = sqrt.(norm.(fv) .* ds)
    restoration_factor = sqrt.(ds ./ norm.(fv))

    for (i,s1) in enumerate(arclengths)
        for (j,s2) in enumerate(arclengths)
            ξ1 = (s1 + s2) / sqrt(2)
            ξ2 = (s1 - s2) / sqrt(2)

            ξ2 = mod(ξ2 + perimeter * sqrt(2) / 4, perimeter * sqrt(2) / 2) - perimeter * sqrt(2) / 4
            ξ1 = mod(ξ1, perimeter * sqrt(2))

            Γ_s[i,j] = sym_factor[i] * sym_factor[j] * etp.(ξ2, ξ1)
        end
        # Enforce particle conservation 
        Γ_s[i,i] -= dot(Γ_s[i, :], restoration_factor) / restoration_factor[i]
    end
    
    open(outfiles[2], "w") do file
        writedlm(file, Γ_s, ",")
    end

end

##### Constants #####
const prec::Float64   = 0.001 # Minimum value of F = f (1 - f) where f is the FD distribution to determine bounds of discretization
const tolerance::Float64  = 0.00001 
const max_iterations::Int = 10000

input_handling()

const band        = ARGS[1]
include(joinpath(@__DIR__, "params", "$(band).jl"))

const temperature = parse(Float64, ARGS[2]) # In units of T/T_F
const num_bins    = parse(Int,     ARGS[3]) # Number of s_values in integration mesh
const perp_num    = parse(Int,     ARGS[4]) # Number of t_values (along v_f) in integration mesh
const row_dim     = 4 * (parse(Int,ARGS[6]) ÷ 4) # Dimension of interpolated output matrix enforced to be divisible by 4

main()