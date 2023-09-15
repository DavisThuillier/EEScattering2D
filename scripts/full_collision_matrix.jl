using EEScattering2D

using DelimitedFiles
using Plots
using LinearAlgebra
using DataFrames
using CSV
import StaticArrays: SVector
using Interpolations

function symmetry_metric(mat::Union{Matrix{Float64}, Matrix{ComplexF64}})
    sym_mat_norm = norm((mat + mat')/2)
    asym_mat_norm = norm((mat - mat')/2)
    symmetry = (sym_mat_norm - asym_mat_norm) / (sym_mat_norm + asym_mat_norm) 

    return 0.5 * (symmetry + 1.0)
end

function matrix_interpolation(fs::Vector{SVector{2,Float64}}, mat::Matrix{Float64}, resized_dim::Int)
    coords = FermiSurfaceMesh.get_arclengths(fs)
    ds = mean(FermiSurfaceMesh.get_ds(fs))

    interpolation_coords = collect(Base._linspace(coords[begin], coords[end], resized_dim + 1))
    pop!(interpolation_coords)
    interpolation_fs = Vector{SVector{2,Float64}}(undef, length(interpolation_coords))
    for i in eachindex(interpolation_coords)
        interpolation_fs[i] = FermiSurfaceMesh.get_momentum(fs, coords, interpolation_coords[i])
    end
    interpolation_ds = FermiSurfaceMesh.get_ds(interpolation_fs)

    interpolation_matrix = Matrix{Float64}(undef, (resized_dim, resized_dim))

    mat = hcat(mat, mat[:, 1])
    mat = vcat(mat, mat[1, :]')
    
    nodes = (coords, coords)
    itp = interpolate(nodes, mat, Gridded(Linear()))

    
    for i in eachindex(interpolation_coords)
        interpolation_matrix[:, i] = sqrt.(interpolation_ds) .* itp.(interpolation_coords, interpolation_coords[i]) * sqrt(interpolation_ds[i]) / ds
        interpolation_fs[i] = FermiSurfaceMesh.get_momentum(fs, coords, interpolation_coords[i])
    end

    return interpolation_matrix, interpolation_fs
end

function main()
    include("params/data_dir.jl")

    fs_filename::String  = joinpath(data_dir, "fermi_surface_$(matrix_dim).csv")
    fermi = CSV.read(fs_filename, DataFrames.DataFrame)
    fs = SVector{2}.(fermi.kx, fermi.ky)


    n_stem = joinpath(data_dir, "Γn_$(matrix_dim)_$(temperature)")
    u_stem = joinpath(data_dir, "Γu_$(matrix_dim)_$(temperature)")
    full_stem = joinpath(data_dir, "Γ_$(matrix_dim)_$(temperature)")

    n_files = String[]
    u_files = String[]
    full_files = String[]
    for file in readdir(data_dir; join=true)
        startswith(file, n_stem) && push!(n_files, file)
        startswith(file, u_stem) && push!(u_files, file)
        startswith(file, full_stem) && push!(full_files, file)
    end
    
    if (length(n_files) != length(u_files) ) 
        println("Error: Missing data")
        exit()
    end

    full_matrix = Matrix{Float64}(undef, matrix_dim, matrix_dim)
    n_matrix = Matrix{Float64}(undef, matrix_dim, matrix_dim)
    u_matrix = Matrix{Float64}(undef, matrix_dim, matrix_dim)
    for i in eachindex(n_files)
        start_index = parse(Int, split(basename(n_files[i]), "_")[4])
        end_index   = parse(Int, split(basename(n_files[i])[begin:(end-4)], "_")[6])

        n_matrix[start_index:end_index, :] = readdlm(n_files[i], ',', Float64)
        u_matrix[start_index:end_index, :] = readdlm(u_files[i], ',', Float64)
    end

    ### Generate Full Matrix ###
    n_matrix[2 + div(matrix_dim, 8): div(matrix_dim, 4) + 1, :] = circshift( reverse( n_matrix[1:div(matrix_dim, 8), :]), (0, 2 * div(matrix_dim, 8) + 1) )
    u_matrix[2 + div(matrix_dim, 8): div(matrix_dim, 4) + 1, :] = circshift( reverse( u_matrix[1:div(matrix_dim, 8), :]), (0, 2 * div(matrix_dim, 8) + 1) )

    n_matrix[1 + div(matrix_dim, 4) : div(matrix_dim, 2), :] = circshift( n_matrix[1 : div(matrix_dim, 4), :], (0, div(matrix_dim, 4)))
    u_matrix[1 + div(matrix_dim, 4) : div(matrix_dim, 2), :] = circshift( u_matrix[1 : div(matrix_dim, 4), :], (0, div(matrix_dim, 4)))

    n_matrix[1 + div(matrix_dim, 2) : matrix_dim, :] = circshift( n_matrix[1 : div(matrix_dim, 2), :], (0, div(matrix_dim, 2)))
    u_matrix[1 + div(matrix_dim, 2) : matrix_dim, :] = circshift( u_matrix[1 : div(matrix_dim, 2), :], (0, div(matrix_dim, 2)))

    full_matrix = ( n_matrix + u_matrix) 

    interpolated_matrix, interpolation_fs = matrix_interpolation(fs, full_matrix, interpolation_dim)

    interpolation_fv = Vector{SVector{2,Float64}}(undef, length(interpolation_fs))
    FermiSurfaceMesh.fill_fermi_velocity!(interpolation_fv, interpolation_fs, hamiltonian)

    @show symmetry_metric(full_matrix)
    @show symmetry_metric(interpolated_matrix)

    # full_matrix = (full_matrix' + full_matrix) / 2 # Symmetrizing matrix
    for i in eachindex(interpolated_matrix[:, 1])
        interpolated_matrix[i,i] -= sum(interpolated_matrix[:, i])# - full_matrix[i,i]
    end

    outfile = "Γ_full_$(interpolation_dim)_$(temperature).csv"
    open(joinpath(data_dir, outfile), "w") do file
        writedlm(file, interpolated_matrix, ",")
    end


    interpolation_fs_file = joinpath(data_dir, "fermi_surface_$(interpolation_dim).csv")
    open(interpolation_fs_file, "w") do file
        println(file, "kx,ky,vx,vy")
        writedlm(file, hcat(first.(interpolation_fs), last.(interpolation_fs), first.(interpolation_fv), last.(interpolation_fv)), ",")
    end
end

include(joinpath(@__DIR__, "params", "gamma.jl"))

main()