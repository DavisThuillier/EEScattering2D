using EEScattering2D

using DelimitedFiles
using Plots
using LinearAlgebra
using DataFrames
using CSV
import StaticArrays: SVector

function main()
    band = "gamma"
    resolution = "100_20"
    matrix_dim = 800
    temperature = 0.01
    run = "2023.07.17"
    
    data_dir = joinpath(@__DIR__, "..", "data", band, resolution, run)

    n_stem = joinpath(data_dir, "Γn_$(matrix_dim)_$(temperature)")
    u_stem = joinpath(data_dir, "Γu_$(matrix_dim)_$(temperature)")

    n_files = String[]
    u_files = String[]
    for file in readdir(data_dir; join=true)
        startswith(file, n_stem) && push!(n_files, file)
        startswith(file, u_stem) && push!(u_files, file)
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

    # Enforce particle conservation by ensuring row sum is null
    for i in 1:(div(matrix_dim, 8) + 1)
        n_matrix[i,i] -= sum(n_matrix[i,:])
        u_matrix[i,i] -= sum(u_matrix[i,:])
    end

    ### Generate Full Matrix ###
    n_matrix[2 + div(matrix_dim, 8): div(matrix_dim, 4) + 1, :] = circshift( reverse( n_matrix[1:div(matrix_dim, 8), :]), (0, 2 * div(matrix_dim, 8) + 1) )
    u_matrix[2 + div(matrix_dim, 8): div(matrix_dim, 4) + 1, :] = circshift( reverse( u_matrix[1:div(matrix_dim, 8), :]), (0, 2 * div(matrix_dim, 8) + 1) )

    n_matrix[1 + div(matrix_dim, 4) : div(matrix_dim, 2), :] = circshift( n_matrix[1 : div(matrix_dim, 4), :], (0, div(matrix_dim, 4)))
    u_matrix[1 + div(matrix_dim, 4) : div(matrix_dim, 2), :] = circshift( u_matrix[1 : div(matrix_dim, 4), :], (0, div(matrix_dim, 4)))

    n_matrix[1 + div(matrix_dim, 2) : matrix_dim, :] = circshift( n_matrix[1 : div(matrix_dim, 2), :], (0, div(matrix_dim, 2)))
    u_matrix[1 + div(matrix_dim, 2) : matrix_dim, :] = circshift( u_matrix[1 : div(matrix_dim, 2), :], (0, div(matrix_dim, 2)))

    u_matrix *= (2pi/ matrix_dim)
    n_matrix *= (2pi/ matrix_dim)
    full_matrix = ( n_matrix + u_matrix) # Prefactor corresponds to the differential in the angular integral when taking the product of full_matrix and a vector
    
    outfile = "Γ_full_$(matrix_dim)_$(temperature).csv"
    open(joinpath(data_dir, outfile), "w") do file
        writedlm(file, full_matrix, ",")
    end
end

main()