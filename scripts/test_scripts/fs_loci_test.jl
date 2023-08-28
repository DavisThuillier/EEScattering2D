using EEScattering2D

import EEScattering2D.FermiSurfaceMesh  as FSM
using Plots
import StaticArrays: SVector


function main(row_dim::Int)
    
    fermi_surface = FSM.generate_fermi_surface(hamiltonian, row_dim)
    perimeter = last(FSM.get_arclengths(fermi_surface))
    loci = [0.0, 0.7, 0.1, 0.8] * perimeter

    colin_width::Float64 = FSM.collinear_width(temperature)
    arclengths = FSM.get_arclengths(fermi_surface)

    new_fs = FSM.get_gaussian_mesh(fermi_surface, loci, num_bins, colin_width / (2 * sqrt(2 * log(2))) )
    new_fv = Vector{SVector{2, Float64}}(undef, length(new_fs))
    FSM.fill_fermi_velocity!(new_fv, new_fs, hamiltonian)

    @show length(new_fs)

    plt = plot(first.(fermi_surface), last.(fermi_surface), aspectratio = 1.0) 
    plot!(plt, first.(new_fs), last.(new_fs), seriestype = :scatter, markershape = :cross, marksize = 0.4, color = :black)


end

include(joinpath("..", "params", "gamma.jl"))
include(joinpath("..", "params", "params.jl"))

main(1000)