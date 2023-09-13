using EEScattering2D

using DelimitedFiles
using Plots
using LinearAlgebra
using DataFrames
using CSV
using LaTeXStrings
import StaticArrays: SVector
import Statistics: mean

function main()
    include("params/data_dir.jl")
    @show data_dir

    start_index::Int = 1
    end_index::Int = 1

    fermi = CSV.read(joinpath(data_dir,"fermi_surface_$(matrix_dim).csv"), DataFrames.DataFrame)
    fs = SVector{2}.(fermi.kx, fermi.ky)
    arclengths = FermiSurfaceMesh.get_arclengths(fs)
    perimeter = last(arclengths)

    @show perimeter - arclengths[end - 1]
    @show arclengths[2]

    filen::String = joinpath(data_dir, "Γn_$(matrix_dim)_$(temperature)_$(start_index)_to_$(end_index).csv")
    fileu::String = joinpath(data_dir, "Γu_$(matrix_dim)_$(temperature)_$(start_index)_to_$(end_index).csv")

    distribution = vec(readdlm(filen, ',', Float64) + readdlm(fileu, ',', Float64))

    distribution[1] -= sum(distribution)

    @show length(distribution)

    plt = plot(arclengths[2:end-1] * 2pi / perimeter, 1e-5 .+ distribution[2:end], proj = :polar)
    display(plt)
end

main()
