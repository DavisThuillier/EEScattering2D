using EEScattering2D

import StaticArrays: SVector
using Plots
using Interpolations
import LinearAlgebra: norm

function main()
    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim)
    momenta, dVs, variance, arclengths = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, 0 + 1, hamiltonian, temperature, prec)

    plt = plot(first.(momenta), last.(momenta), aspectratio = 1.0, seriestype = :scatter, markershape= :cross, markersize = 0.1, color = :green, legend = false)
    display(plt)
end

include("params/params.jl")
include(joinpath(@__DIR__, "params", "$(band).jl"))


main()