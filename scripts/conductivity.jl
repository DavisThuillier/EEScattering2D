using EEScattering2D

using DelimitedFiles
using Plots
using LinearAlgebra
using DataFrames
using CSV
using LaTeXStrings
import StaticArrays: SVector
using Statistics


function main()
    fs_filename::String  = joinpath(data_dir, "fermi_surface_$(matrix_dim).csv")

    if isfile(mat_filename)
        full_matrix::Matrix{Float64} = readdlm(mat_filename, ',', Float64)
        fermi = CSV.read(fs_filename, DataFrames.DataFrame)
        
        fs = SVector{2}.(fermi.kx, fermi.ky)
        ds = FermiSurfaceMesh.get_ds(fs)

        σ = Matrix{Float64}(undef, 2, 2) # Conductivity tensor

        lambdas = reverse(eigvals(full_matrix)) * mean(ds)

        eigenvecs = reverse(eigvecs(full_matrix), dims = 2)

        vx = Vector{ComplexF64}(undef, length(lambdas))
        vy = Vector{ComplexF64}(undef, length(lambdas))
        for i in eachindex(vx)
            # vx[i] = inner_product(fermi.vx, eigenvecs[:, i], fs, hamiltonian, temperature) / (inner_product(eigenvecs[:, i], eigenvecs[:, i], fs, hamiltonian, temperature))
            # vy[i] = inner_product(fermi.vy, eigenvecs[:, i], fs, hamiltonian, temperature) / (inner_product(eigenvecs[:, i], eigenvecs[:, i], fs, hamiltonian, temperature))
            vx[i] = dot(fermi.vx, eigenvecs[:, i]) / dot(eigenvecs[:, i], eigenvecs[:, i])
            vy[i] = dot(fermi.vy, eigenvecs[:, i]) / dot(eigenvecs[:, i], eigenvecs[:, i]) 
        end

        @show abs(vx[3]) / abs(vx[5])

        # Enforce that the overlap with the m = 0 mode is null
        vx[1] = 0.0
        vy[1] = 0.0

        inverse_times = diagm(1 ./ lambdas)
            
        σ[1,1] = real.(inner_product(vx, inverse_times * vx, fs, hamiltonian, temperature))
        σ[1,2] = real.(inner_product(vx, inverse_times * vy, fs, hamiltonian, temperature))
        σ[2,1] = real.(inner_product(vy, inverse_times * vx, fs, hamiltonian, temperature))
        σ[2,2] = real.(inner_product(vy, inverse_times * vy, fs, hamiltonian, temperature))

        prefactor = 193.468 * 4

        println("T = ", temperature)
        println("σ = ", σ)
        println("ρxx = ", 1 / (- (σ[1,1] + σ[2,2])/2 * prefactor))


    else
        println("Data file of full collision matrix does not exist.")
    end
end

include("params/data_dir.jl")
include("params/$(band).jl")

main()