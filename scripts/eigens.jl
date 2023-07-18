using EEScattering2D

using DelimitedFiles
using Plots
using LinearAlgebra
using DataFrames
using CSV
using LaTeXStrings
import StaticArrays: SVector
import Statistics: mean

function mirror_symmetrize(basis::Union{Vector{Vector{ComplexF64}}, Vector{Vector{Float64}}}, mirror_index::Int)
    mirror_mat = Matrix{ComplexF64}(undef, 2, 2)
    for i in 1:2
        for j in 1:2
            mirror_mat[i,j] = dot(basis[i], circshift(reverse(basis[j]), 2 * mirror_index - 1)) 
        end
    end
    weights = eigvecs(mirror_mat)

    w1 = weights[1,1] * basis[1] + weights[2,1] * basis[2]
    w2 = weights[1,2] * basis[1] + weights[2,2] * basis[2]

    theta1 = mean(mod.(angle.(w1), pi))
    theta2 = mean(mod.(angle.(w2), pi))

    w1 = real.(exp(- im * theta1) * w1)
    w2 = real.(exp(- im * theta2) * w2)

    return w1 / norm(w1), w2 / norm(w2) 
end

function fft(v::Union{Vector{ComplexF64}, Vector{Float64}}, m::Int)
    int::ComplexF64 = 0.0
    step::Float64 = 2pi / length(v)

    for k in eachindex(v)
        int += exp( - im * (m * step * (k - 1)) ) * step * v[k]
    end
    return int
end

function main()
    band::String = "gamma"
    resolution::String = "100_20"
    matrix_dim::Int = 800
    temperature::Float64 = 0.01
    run::String = "2023.07.17(2)"

    data_dir = joinpath(@__DIR__, "..", "data", band, resolution, run)
    filename::String = joinpath(data_dir, "Γ_full_$(matrix_dim)_$(temperature).csv")

    if isfile(filename)
        full_matrix::Matrix{Float64} = readdlm(filename, ',', Float64)
        fermi = CSV.read(joinpath(data_dir,"fermi_surface_$(matrix_dim).csv"), DataFrames.DataFrame)
        @show norm(full_matrix' - full_matrix)
        
        fs = SVector{2}.(fermi.kx, fermi.ky)
        fs_norms = norm.(fs)
        thetas = map(x -> mod2pi(atan(x[2], x[1])), fs)
        fv = SVector{2}.(fermi.vx, fermi.vy)

        lambdas = reverse(eigvals(full_matrix))
        eigenvecs = reverse(eigvecs(full_matrix), dims = 2) # Order eigenvectors from smallest to largest

        n = 40
        display(plot(0:n, real.(lambdas[1:n]), color = :green, seriestype = :scatter))

    else
        println("Data file of full collision matrix does not exist.")
        return nothing
    end
    
    scale = 1.0

    println("### Mode Analysis ###")
    maximal_contribution::Float64 = 0.0
    maxima_index::Int = 0
    contribution::Float64 = 0.0

    i = 1
    while i < 20
        println("Eigenvector ", i)
        if abs(real(lambdas[i + 1]) - real(lambdas[i])) < 1e-6
            w1, w2 = mirror_symmetrize([eigenvecs[:, i], eigenvecs[:, i + 1]], div(matrix_dim, 4) + 1)
            i += 1
        else
            println("Single")
            w1 = eigenvecs[:, i] / norm(eigenvecs[:, i])
            w2 = w1
        end

        # normal_contribution = dot(w1, n_matrix * w1) / dot(w1, full_matrix * w1)
        # umklapp_contribution = dot(w1, u_matrix * w1) / dot(w1, full_matrix * w1)
    
        # println("n = ", i/2)
        # @show normal_contribution
        # @show umklapp_contribution

        maximal_contribution = 0.0
        maxima_index = 0
        for j in 0:20
            contribution = abs(fft(w1, j))
            println("m = ", j, ": ", round(contribution, digits = 5))
            if contribution > maximal_contribution
                maximal_contribution = contribution
                maxima_index = j
            end
        end
        println("lambda = ", round(lambdas[i], digits = 5), "; maximal contribution: m = ", maxima_index)
        println()


        plt = plot(thetas, fs_norms .+ scale * real.(w1), title = latexstring("\$ \\lambda = $(round(lambdas[i], digits = 5)) \$"))
        plot!(plt, thetas, fs_norms .+ scale * real(w2), color = :black)
        plot!(plt, thetas, fs_norms, color = :green)
        display(plt)
        i += 1
    end
end

main()