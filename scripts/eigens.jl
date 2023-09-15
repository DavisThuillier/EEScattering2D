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

function display_modes(lambdas, eigenvecs, fs::Vector{SVector{2,Float64}}, fv::Vector{SVector{2,Float64}})
    fs_norms = norm.(fs)
    sqrtspeeds = sqrt.(norm.(fv))

    thetas = map(x -> mod2pi(atan(x[2], x[1])), fs)

    println("### Mode Analysis ###")
    maximal_contribution::Float64 = 0.0
    maximum_index::Int = 0
    contribution::Float64 = 0.0

    even_modes = Vector{Tuple{Int,Int}}(undef, 0)
    odd_modes  = Vector{Tuple{Int,Int}}(undef, 0)

    i = 1
    n = 100
    scale = 1.0
    maximum_index = 0

    while maximum_index != 2
        
        if abs(real(lambdas[i + 1]) - real(lambdas[i])) < 1e-6
            w1, w2 = mirror_symmetrize([eigenvecs[:, i], eigenvecs[:, i + 1]], div(matrix_dim, 4) + 1)
            i += 1
            skip = true
        else
            w1 = eigenvecs[:, i] / norm(eigenvecs[:, i])
            w2 = w1
            skip = false
        end

        maximal_contribution = 0.0
        maximum_index = 0
        for j in 0:500
            contribution = abs(fft(w1, j))
            # println("m = ", j, ": ", round(contribution, digits = 5))
            if contribution > maximal_contribution
                maximal_contribution = contribution
                maximum_index = j
            end
        end

        
        isodd(maximum_index) ? push!(odd_modes, (i, maximum_index)) : push!(even_modes, (i, maximum_index)) 
        if skip
            isodd(maximum_index) ? push!(odd_modes, (i-1, maximum_index)) : push!(even_modes, (i-1, maximum_index))
        end


        println("Eigenvector ", i)
        println("lambda = ", round(lambdas[i] - lambdas[1], digits = 5), "; maximal contribution: m = ", maximum_index)
        if maximum_index in [1, 3, 5]
            @show abs(fft(w1, 1))
            @show abs(fft(w1, 3))
            @show abs(fft(w1, 5))
        end

        println()

        plt = plot(thetas, fs_norms .+ scale * (real.(w1) ./ sqrtspeeds), title = latexstring("\$ \\lambda = $(round(lambdas[i], digits = 5)) , \\mathrm{mode} \\approx $(maximum_index) \$"), proj = :polar)
        #plot!(plt, thetas, fs_norms .+ scale * real(w2), color = :black)
        plot!(plt, thetas, fs_norms, color = :green)
        display(plt)
        i += 1
    end

    
    spectrum = plot(first.(odd_modes), map(x -> lambdas[x], first.(odd_modes)), seriestype = :scatter, label = "Odd Modes")
    plot!(spectrum, first.(even_modes), map(x -> lambdas[x], first.(even_modes)), seriestype = :scatter, label = "Even Modes")
    display(spectrum)
end

function main()
    include("params/data_dir.jl")
    @show data_dir

    if isfile(mat_filename)
        full_matrix::Matrix{Float64} = readdlm(mat_filename, ',', Float64)

        fermi = CSV.read(joinpath(data_dir,"fermi_surface_$(interpolation_dim).csv"), DataFrames.DataFrame)
        fs = SVector{2}.(fermi.kx, fermi.ky)
        fv = SVector{2}.(fermi.vx, fermi.vy)

        arclengths = FermiSurfaceMesh.get_arclengths(fs)
        perimeter  = last(arclengths)
        pop!(arclengths)

        # for i in 1:20:size(full_matrix)[2]
        #     plt = plot(arclengths * 2pi / perimeter, 2e-6 .+ full_matrix[i, :], proj = :polar)
        #     display(plt)
        # end
        # return nothing

        lambdas = real.(reverse(eigvals(full_matrix))) * 5.25e3 # Rates in ps^-1
        eigenvecs = reverse(eigvecs(full_matrix), dims = 2) # Order eigenvectors from smallest to largest

        restoration_factor = sqrt.(norm.(fv) ./ FermiSurfaceMesh.get_ds(fs))
        for i in 1:size(eigenvecs)[2]
            eigenvecs[:, i] = eigenvecs[:, i] .* restoration_factor
        end

        display_modes(lambdas, eigenvecs, fs, fv)


    else
        println("Data file of full collision matrix does not exist.")
        return nothing
    end

    # spectrum2 = plot(last.(odd_modes), -map(x -> lambdas[x], first.(odd_modes)), seriestype = :scatter, label = "Odd Modes", xlims = (0,50))
    # plot!(spectrum2, last.(even_modes), -map(x -> lambdas[x], first.(even_modes)), seriestype = :scatter, label = "Even Modes")
    # plot!(spectrum2, title = "T = $(round(temperature, digits = 4))", ylabel = L"- \lambda_m ", xlabel = "m")
    # display(spectrum2)
end

main()