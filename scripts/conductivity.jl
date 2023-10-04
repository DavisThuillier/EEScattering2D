using EEScattering2D

using DelimitedFiles
using LinearAlgebra
using DataFrames
using CSV
import StaticArrays: SVector
using Statistics


function main()
    fs_filename::String  = joinpath(data_dir, "fermi_surface_$(matrix_dim).csv")

    if isfile(s_matrix_file)
        full_matrix::Matrix{Float64} = readdlm(s_matrix_file, ',', Float64)
        fermi = CSV.read(fs_filename, DataFrames.DataFrame)
        
        fs = SVector{2}.(fermi.kx, fermi.ky)
        ds = mean(FermiSurfaceMesh.get_ds(fs))

        fv = SVector{2}.(fermi.vx, fermi.vy)

        σ = Matrix{ComplexF64}(undef, 2, 2) # Conductivity tensor

        lambdas = reverse(eigvals(full_matrix))
        lambdas = lambdas .- lambdas[1]

        eigenvecs = reverse(eigvecs(full_matrix), dims = 2)

        vx = Vector{ComplexF64}(undef, length(lambdas))
        vy = Vector{ComplexF64}(undef, length(lambdas))
        sqrt_speeds = sqrt.(norm.(fv) ./ ds)
        for i in eachindex(vx)
            vx[i] = dot(fermi.vx ./ sqrt_speeds, eigenvecs[:, i]) / dot(eigenvecs[:, i], eigenvecs[:, i])
            vy[i] = dot(fermi.vy ./ sqrt_speeds, eigenvecs[:, i]) / dot(eigenvecs[:, i], eigenvecs[:, i]) 
        end

        # @show abs.(vx[1:10])
        # @show lambdas[1:5]

        # Enforce that the overlap with the m = 0 mode is null
        vx[1] = 0.0
        vy[1] = 0.0

        inverse_times = diagm(1 ./ lambdas)
        inverse_times[1,1] = 0.0
        # @show findmax(abs.(vx))

        σ[1,1] = real(dot(vx, inverse_times * vx))
        σ[1,2] = real(dot(vx, inverse_times * vy))
        σ[2,1] = real(dot(vy, inverse_times * vx))
        σ[2,2] = real(dot(vy, inverse_times * vy))

        spin_species = 2
        σ = spin_species * σ * 193.468

        # prefactor = 193.468 * 4

        T_F = 6326.35
        println("T = ", temperature * T_F)
        println("σ = ", σ)
        println("ρxx = ", 1 / (- (σ[1,1] + σ[2,2])/2) * 1e11)


    else
        println("Data file of full collision matrix does not exist.")
    end
end

include("params/data_dir.jl")
include("params/$(band).jl")

main()