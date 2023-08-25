using EEScattering2D

import StaticArrays: SVector
using Plots
using Interpolations
import LinearAlgebra: norm
import Statistics: mean, median
using LaTeXStrings
using GeometryBasics
using VoronoiCells

fd(E::Float64, T::Float64) = 1 / (exp(E/T) + 1)

df(E::Float64, T::Float64) = fd(E, T) * (1 - fd(E,T)) / T

function main()
    band == "free" ? (hasBZ = false) : (hasBZ = true)
    temperatures = [0.002, 0.004]
    for temperature in temperatures
        fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim, bz = hasBZ)
        fv = Vector{SVector{2, Float64}}(undef, length(fs))
        FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  
        momenta, dVs, variance, arclengths = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, 100, hamiltonian, temperature, prec; bz = hasBZ)
        
        fd_sum = 0.0
        for i in 2:(size(momenta)[1] - 1)
            fd_sum += df(hamiltonian(momenta[i, 1]), temperature) * norm(momenta[i + 1, 1] - momenta[i-1, 1]) / 2
        end
        @show sqrt(variance) / temperature
        @show fd_sum / norm(fv[1])

        points = map(x -> Point(x[1], x[2]), vec(momenta'))

        BrillouinZone = Rectangle(Point(-0.5, -0.5), Point(0.5, 0.5))
        tess = voronoicells(points, BrillouinZone)
        areas = voronoiarea(tess)

        modulus = size(momenta)[2]
        reduced_cells = tess.Cells[modulus+1:end-modulus]
        #reduced_cells = tess.Cells[1:modulus]
        reduced_points = points[modulus+1:end-modulus]
        #reduced_points = points[1:modulus]
        reduced_tess = Tessellation(reduced_points, BrillouinZone, reduced_cells)
        reduced_areas = voronoiarea(reduced_tess)
        @show sum(reduced_areas)
        @show sum(dVs)

        grid_vals = Base._linspace(-0.5, 0.5, 200)
        grid = collect(Iterators.product(grid_vals, grid_vals))

        plt = heatmap(grid_vals, grid_vals, hamiltonian.(SVector{2}.(grid)), colormap = :nuuk, colorbar_title = latexstring("\$\\epsilon(\\vec{k}) / \\epsilon_F\$"))
        plot!(plt, first.(fs), last.(fs))
        plot!(plt, first.(momenta), last.(momenta), aspectratio = 1.0, seriestype = :scatter, markershape= :cross, markersize = 0.2, color = :black, legend = false, xlims = (-0.5, 0.5), ylims = (-0.5, 0.5))
        plot!(plt, xlabel = L"k_x a_x / 2\pi", ylabel = L"k_y a_y / 2\pi", title = "T = $(temperature)")
        plot!(plt, reduced_tess, fillcolor=:green, linewidth = 0.1)
        display(plt)

        
        
    end
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))


main()