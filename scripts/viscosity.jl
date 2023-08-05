# viscosity.jl
# author: Davis Thuillier
# created: 18 July 2023

using EEScattering2D

include(joinpath(@__DIR__, "params", "gamma.jl"))

function main()
    include("params/data_dir.jl")

    mat_filename::String = joinpath(data_dir, "Γ_full_$(matrix_dim)_$(temperature).csv")
    fs_filename::String  = joinpath(data_dir, "fermi_surface_$(matrix_dim).csv")

    full_matrix::Matrix{Float64} = readdlm(mat_filename, ',', Float64)
    fermi = CSV.read(fs_filename, DataFrames.DataFrame)
        
    fs = SVector{2}.(fermi.kx, fermi.ky)
    ds = FermiSurfaceMesh.get_ds(fs)
    fv = SVector{2}.(fermi.vx, fermi.vy)

    # ηxxxx = dxx.(fs, fv)

end

main()


