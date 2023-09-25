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

function create_files(temperature::Float64, start_index::Int, end_index::Int, num_bins::Int, perp_num::Int)
    data_dir_stem = joinpath(@__DIR__, "..", "data", "$(band)_band", "$(round(temperature, digits = 8))")
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
            println("### Resolution: 192x41")
            # Clears file for writing
        end
    end
    
    return data_dir, filenames
end


function main()
    # data_dir, filenames = create_files(temperature, start_index, end_index, num_bins, perp_num)

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim)
    fv = Vector{SVector{2, Float64}}(undef, length(fs))
    FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  
    # ds = FermiSurfaceMesh.get_ds(fs)
    perimeter = FermiSurfaceMesh.get_perimeter(fs)

    ef::Float64 = 0.55 # Fermi energy in eV
    e0::Float64 = 55.26349406e6 # Vacuum permittivity in e^2 / eV / m
    c::Float64  = 12.68e-10 # z-axis dimension of unit cell in meters

    # Screened Parameters
    alpha = 1 / (ef * e0 * c) # Characteristic non-dimensionalized energy scale for interaction matrix element
    ## prefactor = alpha^2 
    q_squared = alpha * get_dos(fs, fv) / (2pi)^2 # Thomas-Fermi screening wavevector squared 

    # Constant matrix element
    prefactor = (2pi)^2 / get_dos(fs, fv)

    # open(joinpath(data_dir, "fermi_surface_$(length(fs)).csv"), "w") do file
    #     println(file, "kx,ky,vx,vy")
    #     writedlm(file, hcat(first.(fs), last.(fs), first.(fv), last.(fv)), ",")
    # end
    amp_ratio = 8.0

    asym_num = 121
    sym_num = 51

    Γ_ξ = Matrix{Float64}(undef, asym_num, sym_num) 

    asym_mesh = FermiSurfaceMesh.com_gaussian_mesh(-perimeter/2, perimeter/2, asym_num, FermiSurfaceMesh.collinear_width(temperature, perimeter), amp_ratio) / sqrt(2)
    
    # original_coords = Vector{SVector{2,Float64}}(undef, sym_num * length(mesh))
    # k = 1

    sym_mesh = LinRange(0.0,perimeter*sqrt(2), sym_num) # Uniform
    
    s1::Float64 = 0.0
    s2::Float64 = 0.0
    for (j, ξ1) in ProgressBar(enumerate(sym_mesh))
        # plot!(plt, mesh, ξ1 * ones(length(mesh)), seriestype = :scatter, markersize = 0.5, color = :black)
        for (i, ξ2) in enumerate(asym_mesh)
            s1 = mod((ξ1 + ξ2) / sqrt(2), perimeter)
            s2 = mod((ξ1 - ξ2) / sqrt(2), perimeter)
            # @show [s1, s2] / perimeter
            # momenta, _, _, _, _ = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, [s1, s2], hamiltonian, temperature)
            # @show arclengths[loci_indices] / perimeter

            # fig = Figure(fontsize=36, resolution = (1000, 1000))
            # pts = vec(momenta)
            # tri = triangulate(pts)
            # vorn = voronoi(tri)

            # ax = Axis(fig[1,1], title="Full Voronoi Tessellation", xlabel = L"k_x \,(2\pi/a_x)", ylabel = L"k_y \,(2\pi/a_y)")
            # xlims!(ax, -0.5, 0.5)
            # ylims!(ax, -0.5, 0.5)
            # voronoiplot!(ax, vorn, show_generators=true, markercolor=:green, markersize = 4, linewidth = 0.2)

            # resize_to_layout!(fig)
            # display(fig)

            Γ_ξ = prefactor * sum(FermiSurfaceIntegration.contracted_integral(s1, s2, fs, num_bins, perp_num, hamiltonian, temperature, q_squared))
        end
    end

    heatmap(asym_mesh, sym_domain, Γ_ξ)

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
const num_bins    = 120
const perp_num    = 17

main()