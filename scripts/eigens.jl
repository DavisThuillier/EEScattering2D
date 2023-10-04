using EEScattering2D

using DelimitedFiles
using CairoMakie
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
        
        if abs(real(lambdas[i + 1]) - real(lambdas[i])) / abs(real(lambdas[i])) < 1e-2
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

        f = Figure(fontsize = 24)
        ax = PolarAxis(f[1,1], title = latexstring("\$ \\lambda = $(round(lambdas[i], digits = 5)) , \\mathrm{mode} \\approx $(maximum_index) \$"), rlimits = (0.0, 0.55), thetaticklabelsvisible = false)
        lines!(ax, thetas, fs_norms, color = :blue)
        lines!(ax, thetas, fs_norms .+ (real.(w1) ./ sqrtspeeds), color = :green)
        display(f)
        
        i += 1
    end

    spectrum_fig = Figure()
    spectrum_ax = Axis(spectrum_fig[1,1])
    scatter!(spectrum_ax, first.(odd_modes), map(x -> lambdas[x], first.(odd_modes)), label = "Odd Modes")   
    scatter!(spectrum_ax, first.(even_modes), map(x -> lambdas[x], first.(even_modes)), label = "Even Modes")
    axislegend(spectrum_ax)
    display(spectrum_fig) 
    # spectrum = plot(first.(odd_modes), map(x -> lambdas[x], first.(odd_modes)), seriestype = :scatter, label = "Odd Modes")
    # plot!(spectrum, first.(even_modes), map(x -> lambdas[x], first.(even_modes)), seriestype = :scatter, label = "Even Modes")
    # display(spectrum)
end

function main()
    include("params/data_dir.jl")
    @show data_dir

    if isfile(s_matrix_file)
        full_matrix::Matrix{Float64} = readdlm(s_matrix_file, ',', Float64)

        fermi = CSV.read(joinpath(data_dir,"fermi_surface_$(matrix_dim).csv"), DataFrames.DataFrame)
        fs = SVector{2}.(fermi.kx, fermi.ky)
        fv = SVector{2}.(fermi.vx, fermi.vy)

        arclengths = FermiSurfaceMesh.get_arclengths(fs)
        perimeter  = last(arclengths)
        pop!(arclengths)

        fig = Figure(fontsize=36, resolution = (1000, 1000))
        ax = Axis(fig[1,1], xlabel = L"s_1", ylabel = L"s_2", title = latexstring("\$T/T_F = $(round(temperature, digits = 8))\$"))

        scale = ReversibleScale(x -> atan(x), x -> tan(x))
        hm = heatmap!(ax, arclengths, arclengths, full_matrix, colormap = :lisbon, colorscale = scale, colorrange = (-0.5e-6,1e-6))
        Colorbar(fig[:,end+1], hm)
        
        display(fig)

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
    
end

main()