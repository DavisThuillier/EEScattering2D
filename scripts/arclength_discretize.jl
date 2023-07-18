using EEScattering2D

import StaticArrays: SVector
using Plots
using Interpolations
import LinearAlgebra: norm
import Statistics: mean

function main()
    band == "free" ? (hasBZ = false) : (hasBZ = true)
    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim, bz = hasBZ)
    momenta, dVs, variance, arclengths = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, 1, hamiltonian, temperature, prec; bz = hasBZ)

    plt = plot(first.(momenta), last.(momenta), aspectratio = 1.0, seriestype = :scatter, markershape= :cross, markersize = 0.1, color = :black, legend = false)
    display(plt)
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))


main()