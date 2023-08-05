using EEScattering2D

using DelimitedFiles
using Plots
using LinearAlgebra
using DataFrames
using CSV
using LaTeXStrings
import StaticArrays: SVector
using Statistics

# function inner_product(v1::Union{Vector{Float64}, Vector{ComplexF64}}, v2::Union{Vector{Float64}, Vector{ComplexF64}}, weights::Vector{Float64})
#     length(v1) == length(v2) || throw(DimensionMismatch("vectors in inner product must have same dimensions"))
#     product = 0.0
#     for i in eachindex(v1)
#         product += weights[i] * conj(v1[i]) * v2[i]
#     end

#     return product
# end

function main()
    
    fs_filename::String  = joinpath(data_dir, "fermi_surface_$(matrix_dim).csv")

    if isfile(mat_filename)
        full_matrix::Matrix{Float64} = readdlm(mat_filename, ',', Float64)
        fermi = CSV.read(fs_filename, DataFrames.DataFrame)
        
        fs = SVector{2}.(fermi.kx, fermi.ky)
        ds = FermiSurfaceMesh.get_ds(fs)

        fv = SVector{2}.(fermi.vx, fermi.vy)

        σ = Matrix{Float64}(undef, 2, 2) # Conductivity tensor

        lambdas = reverse(eigvals(full_matrix))
        eigenvecs = reverse(eigvecs(full_matrix), dims = 2)

        vx = Vector{ComplexF64}(undef, length(lambdas))
        vy = Vector{ComplexF64}(undef, length(lambdas))
        for i in eachindex(vx)
            # vx[i] = inner_product(fermi.vx, eigenvecs[:, i], ds) / (inner_product(eigenvecs[:, i], eigenvecs[:, i], ds))
            vx[i] = dot(fermi.vx, eigenvecs[:, i]) / dot(eigenvecs[:, i], eigenvecs[:, i])
            vy[i] = dot(fermi.vy, eigenvecs[:, i]) / dot(eigenvecs[:, i], eigenvecs[:, i])
            # vy[i] = inner_product(fermi.vy, eigenvecs[:, i], ds) / (inner_product(eigenvecs[:, i], eigenvecs[:, i], ds))
        end
        
    σ[1,1] = real.(inner_product(vx, diagm(lambdas) * vx, fs, hamiltonian, temperature))
    σ[1,2] =  real.(inner_product(vx, diagm(lambdas) * vy, fs, hamiltonian, temperature))
    σ[2,1] = real.(inner_product(vy, diagm(lambdas) * vx, fs, hamiltonian, temperature))
    σ[2,2] = real.(inner_product(vy, diagm(lambdas) * vy, fs, hamiltonian, temperature))
        display(σ * 5.25e3)

    else
        println("Data file of full collision matrix does not exist.")
    end
end

include("params/data_dir.jl")
include("params/$(band).jl")

main()