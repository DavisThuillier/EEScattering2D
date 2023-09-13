using EEScattering2D

import StaticArrays: SVector
using Plots
import LinearAlgebra: norm
import Statistics: mean, median, max
using LaTeXStrings
using GeometryBasics
using VoronoiCells


fd(E::Float64, T::Float64) = 1 / (exp(E/T) + 1)

function main()
    band == "free" ? (hasBZ = false) : (hasBZ = true)
    Tf = 6326.35178 # K
    temperatures = collect(4:4) / Tf

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim, bz = hasBZ)
    arclengths = FermiSurfaceMesh.get_arclengths(fs)


    fs_plot = plot(first.(fs), last.(fs), aspect_ratio = 1.0)
    display(fs_plot)
    fv = Vector{SVector{2, Float64}}(undef, length(fs))
    FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  

    
    # grid_vals = Base._linspace(-0.5, 0.5, 200)
    # grid = collect(Iterators.product(grid_vals, grid_vals))

    # rect = Rectangle(Point(-0.5, -0.5), Point(0.5, 0.5))

    # for temperature in temperatures
    #     @show temperature
    #     loci = [0.1, 0.8] * perimeter
    #     @time momenta, dVs, variance, arclengths, loci_indices = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, loci, hamiltonian, temperature, prec; bz = hasBZ)

    #     size(unique(momenta)) != size(vec(momenta)) && println("Nonunique!")

    #     t_num = size(momenta)[1]  
    #     s_num = size(momenta)[2]   

    #     points = map(x -> Point(x[1], x[2]), vec(momenta'))
    #     tess = voronoicells(points, rect)
    #     areas = voronoiarea(tess)
    #     areas = reshape(areas, size(momenta))
    #     fill!(areas[begin, :], 0.0)
    #     fill!(areas[end, :], 0.0)

    #     modulus = size(momenta)[2]
    #     reduced_cells = tess.Cells[modulus+1:end-modulus]
    #     reduced_points = points[modulus+1:end-modulus]
        
    #     reduced_tess = Tessellation(reduced_points, rect, reduced_cells)

    #     # plt = plot(first.(momenta), last.(momenta), aspectratio = 1.0, seriestype = :scatter, markershape= :cross, markersize = 0.2, color = :black, legend = false, xlims = (-0.55, 0.55), ylims = (-0.55, 0.55))
    #     # plot!(plt, xlabel = L"k_x a_x / 2\pi", ylabel = L"k_y a_y / 2\pi", title = "T = $(temperature)")
    #     # plot!(plt, reduced_tess, color=:green, linewidth = 0.1)
    #     # display(plt) 
        
    # end
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))


main()