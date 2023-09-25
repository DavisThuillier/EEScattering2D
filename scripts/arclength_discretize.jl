using EEScattering2D

import StaticArrays: SVector
import LinearAlgebra: norm
import Statistics: mean, median, max
using LaTeXStrings
using GeometryBasics
using DelaunayTriangulation, CairoMakie


fd(E::Float64, T::Float64) = 1 / (exp(E/T) + 1)

function main()
    band == "free" ? (hasBZ = false) : (hasBZ = true)
    Tf = 6326.35178 # K
    temperatures = collect(8:4:20) / Tf

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim, bz = hasBZ)
    arclengths = FermiSurfaceMesh.get_arclengths(fs)
    perimeter = last(arclengths)

    num_bins = 192
    perp_num = 13

    fv = Vector{SVector{2, Float64}}(undef, length(fs))
    FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  

    grid_vals = Base._linspace(-0.5, 0.5, 200)
    grid = collect(Iterators.product(grid_vals, grid_vals))

    # rect = Rectangle(Point(-0.5, -0.5), Point(0.5, 0.5))

    for temperature in temperatures
        @show temperature
        loci = [0.0, 0.1] * perimeter
        @time momenta, dVs, variance, arclengths, loci_indices = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, loci, hamiltonian, temperature, prec; bz = hasBZ)
        @show length(arclengths)

        size(unique(momenta)) != size(vec(momenta)) && println("Nonunique!")

        t_num = size(momenta)[1]  
        s_num = size(momenta)[2]   

        fig = Figure(fontsize=24)
        pts = vec(momenta)
        
        # inner_boundary = map(x -> (x[1], x[2]), momenta[begin, :])
        # push!(inner_boundary, inner_boundary[begin])
        # outer_boundary = reverse(map(x -> (x[1], x[2]), momenta[end, :]))
        # push!(outer_boundary, outer_boundary[begin])

        # nodes, pts = convert_boundary_points_to_indices([[outer_boundary]]; existing_points = pts)
        # cons_tri = triangulate(pts; boundary_nodes = nodes, check_arguments = false)
        tri = triangulate(pts)
        vorn = voronoi(tri)
        gens = Vector{SVector{2, Float64}}(undef, length(pts))
        for i in eachindex(gens)
            x = DelaunayTriangulation.get_generator(vorn, i)
            gens[i] = SVector{2}(x[1], x[2])
        end
        @show gens == pts
        # reduced_vorn = VoronoiTessellation(tri, )
        # areas = get_area.((vorn,), 1:length(pts))
        # dVs = reshape(areas, size(momenta))
        # for i in 1:size(dVs)[2]
        #     dVs[begin, i] = 0.0
        #     dVs[end, i]   = 0.0
        # end

        # @show dVs[begin, :]
        ax = Axis(fig[1, 1], title="Clipped Voronoi tessellation", titlealign=:left, width=1000, height=1000)
        voronoiplot!(ax, vorn, show_generators=true, color=:white, markersize = 4, markercolor = :green)

        resize_to_layout!(fig)
        display(fig)


        # points = map(x -> Point(x[1], x[2]), vec(momenta'))
        # tess = voronoicells(points, rect)
        # areas = voronoiarea(tess)
        # areas = reshape(areas, size(momenta))
        # fill!(areas[begin, :], 0.0)
        # fill!(areas[end, :], 0.0)

        # modulus = size(momenta)[2]
        # reduced_cells = tess.Cells[modulus+1:end-modulus]
        # reduced_points = points[modulus+1:end-modulus]
        
        # reduced_tess = Tessellation(reduced_points, rect, reduced_cells)

        # plt = plot(first.(momenta), last.(momenta), aspectratio = 1.0, seriestype = :scatter, markershape= :cross, markersize = 0.2, color = :black, legend = false, xlims = (-0.55, 0.55), ylims = (-0.55, 0.55))
        # plot!(plt, xlabel = L"k_x a_x / 2\pi", ylabel = L"k_y a_y / 2\pi", title = "T = $(temperature)")
        # # plot!(plt, reduced_tess, color=:green, linewidth = 0.1)
        # display(plt) 
        
    end
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))


main()