using EEScattering2D

using DelimitedFiles
using Plots
using LinearAlgebra
using DataFrames
using CSV
using LaTeXStrings
import StaticArrays: SVector
using Statistics

function inner_product(v1::Union{Vector{Float64}, Vector{ComplexF64}}, v2::Union{Vector{Float64}, Vector{ComplexF64}}, weights::Vector{Float64})
    length(v1) == length(v2) || throw(DimensionMismatch("vectors in inner product must have same dimensions"))
    product = 0.0
    for i in eachindex(v1)
        product += weights[i] * conj(v1[i]) * v2[i]
    end

    return product
end

function main()
    band::String = "gamma"
    resolution::String = "100_20"
    matrix_dim::Int = 800
    temperature::Float64 = 0.01
    run::String = "2023.07.17"

    data_dir = joinpath(@__DIR__, "..", "data", band, resolution, run)
    mat_filename::String = joinpath(data_dir, "Γ_full_$(matrix_dim)_$(temperature).csv")
    fs_filename::String  = joinpath(data_dir, "fermi_surface_$(matrix_dim).csv")

    if isfile(mat_filename)
        full_matrix::Matrix{Float64} = readdlm(mat_filename, ',', Float64)
        @show norm(full_matrix' - full_matrix)
        fermi = CSV.read(fs_filename, DataFrames.DataFrame)
        
        fs = SVector{2}.(fermi.kx, fermi.ky)
        ds = FermiSurfaceMesh.get_ds(fs)
        fv = SVector{2}.(fermi.vx, fermi.vy)

        σ = Matrix{Float64}(undef, 2, 2) # Conductivity tensor

        lambdas = reverse(eigvals(full_matrix))
        eigenvecs = reverse(eigvecs(full_matrix), dims = 2)

        n = 8
        inner_products = Matrix{ComplexF64}(undef, n, n)
        for i in 1:n
            for j in 1:n
                # inner_products[i,j] = inner_product(eigenvecs[:, i], eigenvecs[:, j], ds)
                inner_products[i,j] = dot(eigenvecs[:, i], eigenvecs[:, j])
            end
        end
        display(inner_products)
        return nothing


        cutoff = 10
        b = Vector{ComplexF64}(undef, cutoff)
        for i in eachindex(b)
            b[i] = inner_product(fermi.vy, eigenvecs[:, i], ds) / (inner_product(eigenvecs[:, i], eigenvecs[:, i], ds))
        end

        # @show inner_product(fermi.vx, eigenvecs[:, 2], ds) / inner_product(eigenvecs[:, 2], eigenvecs[:, 2], ds)

        vy_approx = zeros(ComplexF64, matrix_dim)
        for i in eachindex(b)
            vy_approx += b[i] * eigenvecs[:, i]
        end
        plt = plot(real.(vy_approx))
        plot!(plt, fermi.vy)
        display(plt)

    else
        println("Data file of full collision matrix does not exist.")
    end
end

main()