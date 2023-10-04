using EEScattering2D

import StaticArrays: SVector
import LinearAlgebra: norm, dot
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

    filenames = map(x -> "Γ_$(x)_$(row_dim)_$(round(temperature, digits = 8)).csv", ["ξ", "s"]) 

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
    data_dir, outfiles = create_files(temperature, num_bins, perp_num)

    # Enforce that the dimension of the matrix is divisible by 4
    dim = 4 * (row_dim ÷ 4)

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, dim) # Uniform Fermi Surface
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

    fig = Figure(fontsize=36, resolution = (1000, 1000))
    ax = Axis(fig[1,1], xlabel = L"ξ_2", ylabel = L"ξ_1", title = latexstring("\$T/T_F = $(round(temperature, digits = 8))\$"))
    hm = heatmap!(ax, Γ_ξ, colormap = Reverse(:davos), colorrange = (-0.002,0.001))
    Colorbar(fig[:,end+1], hm)
    
    display(fig)

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
        Γ_s[i,i] -= dot(Γ_s[i, :], restoration_factor) / restoration_factor[i]
    end

    fig = Figure(fontsize=36, resolution = (1000, 1000))
    ax = Axis(fig[1,1], xlabel = L"s_1", ylabel = L"s_2", title = latexstring("\$T/T_F = $(round(temperature, digits = 8))\$"))
    hm = heatmap!(ax, Γ_s, colormap = Reverse(:davos), colorrange = (-0.002,0.001))
    Colorbar(fig[:,end+1], hm)
    
    display(fig)

    ##################
    ## Symmetrizing ##
    ##################

    # for i in 1:row_dim
    #     for j in 1:row_dim
    #         Γ_s[i,j] = sym_factor[i] * sym_factor[j] * Γ_s[i,j]
    #     end
    #     Γ_s[i,i] = dot(Γ_s[:, i], sqrt.(ds ./ norm.(fv))) * sqrt(norm(fv[i]) / ds[i]) # Enforcing particle conservation
    # end

    
    open(outfiles[2], "w") do file
        writedlm(file, Γ_s, ",")
    end

end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))

# input_handling()

# const temperature = parse(Float64, ARGS[1])
# # const start_index = parse(Int,     ARGS[2])
# # const end_index   = parse(Int,     ARGS[3])
# const num_bins    = parse(Int,     ARGS[2])
# const perp_num    = parse(Int,     ARGS[3])
const temperature::Float64 = 0.001581
const num_bins::Int        = 75
const perp_num::Int        = 17
const row_dim::Int         = 1200

main()