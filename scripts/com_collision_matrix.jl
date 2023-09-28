using EEScattering2D

import StaticArrays: SVector
import LinearAlgebra: norm
using Interpolations
using DelimitedFiles
using ProgressBars

using DelaunayTriangulation, CairoMakie 
using LaTeXStrings

function input_handling()
    if length(ARGS) != 3
        println("Insufficient number of arguments.")
        println("Format: com_collision_matrix.jl [T::F64] [n_s::Int] [n_t::Int]")
        exit()
    end
end

function create_files(temperature::Float64, num_bins::Int, perp_num::Int)
    data_dir_stem = joinpath(@__DIR__, "..", "data", "$(band)_band", "$(round(temperature, digits = 8))", "$(num_bins)_$(perp_num)")
    data_dir = data_dir_stem

    Γ_outfile = joinpath(data_dir, "Γ_full_$(row_dim)_$(round(temperature, digits = 8)).csv")
    
    # Check if any of the output files exist in the desired directory
    j = 0
    while sum(isfile.(joinpath.(data_dir, Γ_outfile))) > 0
        j += 1
        data_dir = data_dir_stem * "($(j))"
        !isdir(data_dir) && break
    end
    !isdir(data_dir) && mkpath(data_dir)

    open(Γ_outfile, "w") do file
        # Clear output file for writing 
    end
    
    return data_dir, Γ_outfile
end

function main()
    data_dir, outfile = create_files(temperature, num_bins, perp_num)

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim) # Uniform Fermi Surface
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

    ###########################
    ## Center-of-mass meshes ##
    ###########################

    sym_num = 2*div(num_bins,4) + 1 # Enforce that sym_num is odd 
    asym_mesh = FermiSurfaceMesh.com_gaussian_mesh(-perimeter/2, perimeter/2, num_bins, FermiSurfaceMesh.collinear_width(temperature, perimeter), amp_ratio) / sqrt(2) # Mesh for ξ2 ≡ (s1 - s2) / √2
    sym_mesh = collect(LinRange(0.0,perimeter*sqrt(2), sym_num))  # Mesh for ξ1 ≡ (s1 + s2) / √2
    asym_num = length(asym_mesh)
    
    s1::Float64 = 0.0
    s2::Float64 = 0.0
    Γ_ξ = Matrix{Float64}(undef, asym_num, sym_num) 
    for (j, ξ1) in ProgressBar(enumerate(sym_mesh))
        for (i, ξ2) in enumerate(asym_mesh)
            s1 = mod((ξ1 + ξ2) / sqrt(2), perimeter)
            s2 = mod((ξ1 - ξ2) / sqrt(2), perimeter)
            
            Γ_ξ[i,j] = prefactor * sum(FermiSurfaceIntegration.contracted_integral(s1, s2, fs, num_bins, perp_num, hamiltonian, temperature, q_squared)) # Sums the normal and umklapp contributions
        end
    end
    ###################
    ## Interpolation ##
    ###################
    itp = interpolate((asym_mesh, sym_mesh), Γ_ξ, Gridded(Linear()))

    # Uniform square mesh for ξ1, ξ2
    # Twice as many points are used for ξ1, since the domain of ξ1
    # is twice as large as the domain of ξ2
    ξ2_mesh = LinRange(-perimeter/2,perimeter/2,row_dim) / sqrt(2)
    ξ1_mesh = LinRange(0.0,perimeter,2*row_dim) * sqrt(2)

    # Collision matrix in (s1,s2) coordinates
    # s1 ↦ row 
    # s2 ↦ column
    Γ_s = Matrix{Float64}(undef, row_dim, row_dim)
    s_mesh = LinRange(0.0, perimeter, row_dim)

    for ξ1 in ξ1_mesh
        for ξ2 in ξ2_mesh
            s1 = mod((ξ1 + ξ2) / sqrt(2), perimeter)
            s2 = mod((ξ1 - ξ2) / sqrt(2), perimeter)

            _, k1 = findmin(abs.(mod.(s1 .- s_mesh,perimeter)))
            _, k2 = findmin(abs.(mod.(s2 .- s_mesh,perimeter)))
            Γ_s[k1, k2] = itp.(ξ2, ξ1)
        end
    end

    ##################
    ## Symmetrizing ##
    ##################

    sym_factor = sqrt.(norm.(fv) .* ds)

    @show size(Γ_s)
    
    open(outfile, "w") do file
        writedlm(file, Γ_s, ",")
    end


    # Γ_s = prefactor * (sym_factor' .* Γ_s * sym_factor)

    # 
    # open(outfile, "w") do file
    #     writedlm(file, Γ_s, ",")
    # end

end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))

# input_handling()

# const temperature = parse(Float64, ARGS[1])
# # const start_index = parse(Int,     ARGS[2])
# # const end_index   = parse(Int,     ARGS[3])
# const num_bins    = parse(Int,     ARGS[2])
# const perp_num    = parse(Int,     ARGS[3])
const temperature = 0.001581
const num_bins    = 50
const perp_num    = 13

main()