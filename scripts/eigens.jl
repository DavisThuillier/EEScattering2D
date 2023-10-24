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
    T_f = 6326.35

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

    while i < 200
        
        if abs(real(lambdas[i + 1]) - real(lambdas[i])) / abs(real(lambdas[i])) < 1e-3
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

        
        if isodd(maximum_index)
            !(maximum_index ∈ last.(odd_modes)) && push!(odd_modes, (i, maximum_index))
        else
            !(maximum_index ∈ last.(even_modes)) && push!(even_modes, (i, maximum_index))
        end
        if skip
            if isodd(maximum_index)
                !(maximum_index ∈ last.(odd_modes)) && push!(odd_modes, (i-1, maximum_index))
            else
                !(maximum_index ∈ last.(even_modes)) && push!(even_modes, (i-1, maximum_index))
            end
        end


        println("Eigenvector ", i)
        println("lambda = ", round(lambdas[i], digits = 5), "; maximal contribution: m = ", maximum_index)
        if isodd(maximum_index) && i < 17
            @show fft(w1, 1)
            @show fft(w1, 3)
            @show fft(w1, 5)
        end
        println()

        # f = Figure(fontsize = 24)
        # ax = PolarAxis(f[1,1], title = latexstring("\$ \\lambda = $(round(-lambdas[i], digits = 5)) \\,\\mathrm{ps}^{-1}, m \\approx $(maximum_index) \$"), rlimits = (0.0, 0.55), thetaticklabelsvisible = false)
        # lines!(ax, thetas, fs_norms, color = :blue)
        # lines!(ax, thetas, fs_norms .+ (real.(w1) ./ sqrtspeeds), color = :green)
        # display(f)
        
        i += 1
    end

    spectrum_fig = Figure(fontsize = 24)
    spectrum_ax = Axis(spectrum_fig[1,1],title = latexstring("\$\\gamma\$ Band, \$(T = $(round(temperature * T_f, digits = 2)) \\, \\mathrm{K})\$"), ylabel = L"\lambda (\mathrm{ps}^{-1})", xticksvisible = false)
    scatter!(spectrum_ax, first.(odd_modes), map(x -> -lambdas[x], first.(odd_modes)), label = "Odd Mode")   
    scatter!(spectrum_ax, first.(even_modes), map(x -> -lambdas[x], first.(even_modes)), label = "Even Mode")
    axislegend(spectrum_ax, valign = :bottom)
    display(spectrum_fig) 

    odd_even_fig = Figure(fontsize = 24)
    odd_even_ax = Axis(odd_even_fig[1,1], ylabel = L"\lambda (\mathrm{ps}^{-1})", xlabel = L"m")
    if band == "gamma"
        odd_even_ax.title = latexstring("\$\\gamma\$ Band, \$(T = $(round(temperature * T_f, digits = 2)) \\, \\mathrm{K})\$")
        odd_even_ax.xlabel = latexstring("Maximal Harmonic Contribution \$(m)\$")
    elseif band == "free"
        title = latexstring("2DEG \$(T = $(temperature) T_F)\$")
    end

    scatter!(odd_even_ax, last.(odd_modes), map(x -> -lambdas[x], first.(odd_modes)), label = "Odd Mode")   
    scatter!(odd_even_ax, last.(even_modes), map(x -> -lambdas[x], first.(even_modes)), label = "Even Mode")
    axislegend(odd_even_ax, valign = :bottom)
    display(odd_even_fig) 
end

function main()
    include("params/data_dir.jl")
    @show data_dir
    T_f = 6326.35

    if isfile(s_matrix_file)
        full_matrix::Matrix{Float64} = readdlm(s_matrix_file, ',', Float64)

        fermi = CSV.read(joinpath(data_dir,"fermi_surface_$(matrix_dim).csv"), DataFrames.DataFrame)
        fs = SVector{2}.(fermi.kx, fermi.ky)
        fv = SVector{2}.(fermi.vx, fermi.vy)

        # n = 9
        # weights = zeros(ComplexF64, n)
        # normv = norm(fermi.vx)
        # for m in 1:n
        #     weights[m] = fft(fermi.vx / normv, m)
        #     @show round(weights[m], digits = 6)
        # end
    
        arclengths = FermiSurfaceMesh.get_arclengths(fs)
        perimeter  = last(arclengths)
        pop!(arclengths)

        fig = Figure(fontsize=30, resolution = (1000, 1000))
        ax = Axis(fig[1,1], xlabel = L"s_1 / (2\pi / a)", ylabel = L"s_2 / (2\pi / a)")
        ax.title = latexstring("\$T = $(round(temperature * T_f, digits = 2)) \\, \\mathrm{K}\$")

        scale = ReversibleScale(x -> atan(x), x -> tan(x))
        hm = heatmap!(ax, arclengths, arclengths, full_matrix * 5.25e3, colormap = :lisbon, colorscale = scale, colorrange = (-3e-4,9e-4))
        Colorbar(fig[:,end+1], hm, label = L"̌\check \xi", labelrotation = pi/2)
        
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