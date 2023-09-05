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
    temperatures = [0.002]
    for temperature in temperatures
        fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim, bz = hasBZ)
        fv = Vector{SVector{2, Float64}}(undef, length(fs))
        FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  

        perimeter = last(FermiSurfaceMesh.get_arclengths(fs))
        grid_vals = Base._linspace(-0.5, 0.5, 200)
        grid = collect(Iterators.product(grid_vals, grid_vals))

        loci = [0.1] * perimeter
        momenta, dVs, variance, arclengths, loci_indices = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, loci, hamiltonian, temperature, prec; bz = hasBZ)
        

        rect = Rectangle(Point2(-0.51,-0.51), Point2(0.51,0.51))
        points = map(x -> Point(x[1], x[2]), vec(unique(momenta)))

        t_num = size(momenta)[1]  
        s_num = size(momenta)[2]
        for i in 1:10:s_num    
            loci = [arclengths[begin], arclengths[i]]    

            mesh_1, _, variance, _ , loci_indices_1 = FermiSurfaceMesh.discretize(fs, s_num, t_num, loci, hamiltonian, temperature, prec)
            sizeof(unique(mesh_1)) != sizeof(mesh_1) && println("Nonunique")

            points = map(x -> Point(x[1], x[2]), vec(mesh_1'))
            tess = voronoicells(points, rect)
            modulus = size(mesh_1)[2]
            reduced_cells = tess.Cells[modulus+1:end-modulus]
            reduced_points = points[modulus+1:end-modulus]
            
            reduced_tess = Tessellation(reduced_points, rect, reduced_cells)

            plt = plot(first.(mesh_1), last.(mesh_1), aspectratio = 1.0, seriestype = :scatter, markershape= :cross, markersize = 0.2, color = :black, legend = false, xlims = (-0.5, 0.5), ylims = (-0.5, 0.5))
            plot!(plt, xlabel = L"k_x a_x / 2\pi", ylabel = L"k_y a_y / 2\pi", title = "T = $(temperature)")
            plot!(plt, reduced_tess, color=:green, linewidth = 0.1)
            display(plt)
        end     
        
    end
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))


main()