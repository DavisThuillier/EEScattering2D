using EEScattering2D

import StaticArrays: SVector
import LinearAlgebra: norm
import Statistics: mean, median, max
using LaTeXStrings
using ProgressBars

using DelaunayTriangulation, CairoMakie


fd(E::Float64, T::Float64) = 1 / (exp(E/T) + 1)

function main()
    band == "free" ? (hasBZ = false) : (hasBZ = true)
    Tf = 6326.35178 # K
    # temperature = collect(8:4:20) / Tf
    temperature = 0.002

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, 1200, bz = hasBZ)
    arclengths = FermiSurfaceMesh.get_arclengths(fs)
    perimeter = last(arclengths)

    num_bins = 208
    perp_num = 41

    # fv = Vector{SVector{2, Float64}}(undef, length(fs))
    # FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  

    grid_vals = Base._linspace(-0.5, 0.5, 200)
    grid = map(x -> SVector{2}(x[1], x[2]), collect(Iterators.product(grid_vals, grid_vals)) )
    energies = hamiltonian.(grid)

    # amp_ratio = 8.0

    # sym_num = 2*div(num_bins,4) + 1 # Enforce that sym_num is odd 

    # asym_mesh = FermiSurfaceMesh.com_gaussian_mesh(-perimeter/2, perimeter/2, num_bins, FermiSurfaceMesh.collinear_width(temperature, perimeter), amp_ratio) / sqrt(2)
    # sym_mesh = collect(LinRange(3 * perimeter / 4,perimeter, sym_num)) * sqrt(2) # Uniform

    # asym_num = length(asym_mesh)
    
    # s1::Float64 = 0.0
    # s2::Float64 = 0.0
    # for (j, ξ1) in ProgressBar(enumerate(sym_mesh))
    #     for (i, ξ2) in enumerate(asym_mesh)
    #         s1 = mod((ξ1 + ξ2) / sqrt(2), perimeter)
    #         s2 = mod((ξ1 - ξ2) / sqrt(2), perimeter)
            
    #         momenta, _, _, arclengths, loci_indices = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, [s1, s2], hamiltonian, temperature)
    #         # @show ([s1, s2] - arclengths[loci_indices], perimeter) / perimeter
    #         println(([s1, s2] - mod.(arclengths[loci_indices], perimeter)) / perimeter)

    #         fig = Figure(fontsize=36, resolution = (1000, 1000))
    #         pts = vec(momenta)
    #         tri = triangulate(pts)
    #         vorn = voronoi(tri)

    #         ax = Axis(fig[1,1], title="Full Voronoi Tessellation", xlabel = L"k_x \,(2\pi/a_x)", ylabel = L"k_y \,(2\pi/a_y)")
    #         heatmap!(ax, grid_vals, grid_vals, energies)
    #         xlims!(ax, -0.5, 0.5)
    #         ylims!(ax, -0.5, 0.5)
    #         voronoiplot!(ax, vorn, show_generators=true, markercolor=:green, markersize = 4, linewidth = 0.2)
            

    #         resize_to_layout!(fig)
    #         display(fig)
            
    #     end
    # end

    momenta, _, _, arclengths, loci_indices = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, [0.0, 0.2] * perimeter, hamiltonian, temperature, bz = hasBZ)
    f = Figure(fontsize=36, resolution = (1000, 1000))
    ax = Axis(f[1,1], xlabel = L"k_x \,(2\pi/a_x)", ylabel = L"k_y \,(2\pi/a_y)")

    heatmap!(ax, grid_vals, grid_vals, energies, colormap = Reverse(:lisbon), colorrange = (-4, 1))
    plot!(ax, first.(fs), last.(fs), linewidth = 0.5, color = :green)
    plot!(ax, vec(first.(momenta)), vec(last.(momenta)), markersize = 4.0)

    display(f)
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))


main()