using EEScattering2D

import StaticArrays: SVector
using Plots
using Interpolations
import LinearAlgebra: norm
import Statistics: mean
using LaTeXStrings

function main()
    band == "free" ? (hasBZ = false) : (hasBZ = true)
    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim, bz = hasBZ)
    momenta, dVs, variance, arclengths = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, 20, hamiltonian, temperature, prec; bz = hasBZ)

    grid_vals = Base._linspace(-0.5, 0.5, 200)
    grid = collect(Iterators.product(grid_vals, grid_vals))

    plt = heatmap(grid_vals, grid_vals, hamiltonian.(SVector{2}.(grid)), colormap = :nuuk, colorbar_title = latexstring("\$\\epsilon(\\vec{k}) / \\epsilon_F\$"))
    plot!(plt, first.(momenta), last.(momenta), aspectratio = 1.0, seriestype = :scatter, markershape= :cross, markersize = 0.2, color = :black, legend = false, xlims = (-0.5, 0.5), ylims = (-0.5, 0.5))
    plot!(plt, xlabel = L"k_x a_x / 2\pi", ylabel = L"k_y a_y / 2\pi")
    display(plt)
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))


main()