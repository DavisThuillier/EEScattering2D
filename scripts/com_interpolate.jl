using EEScattering2D

using DelimitedFiles
using CSV
using DataFrames
import LinearAlgebra: dot, norm
import StaticArrays: SVector
using LaTeXStrings

using CairoMakie

function plot_matrix(Γ, s)
    fig = Figure(fontsize=36, resolution = (1000, 1000))
    ax = Axis(fig[1,1], xlabel = L"s_1", ylabel = L"s_2", title = latexstring("\$T/T_F = $(round(temperature, digits = 8))\$"))
    hm = heatmap!(ax, s, s,Γ, colormap = Reverse(:davos))
    Colorbar(fig[:,end+1], hm)
    
    display(fig)
end

function main()
    include("params/data_dir.jl")

    fs_filename::String  = joinpath(data_dir, "fermi_surface_$(matrix_dim).csv")
    fermi = CSV.read(fs_filename, DataFrames.DataFrame)
    fs = SVector{2}.(fermi.kx, fermi.ky)
    arclengths = FermiSurfaceMesh.get_arclengths(fs)
    pop!(arclengths)
    fv = SVector{2}.(fermi.vx, fermi.vy)
    ds = FermiSurfaceMesh.get_ds(fs)


    Γ::Matrix{Float64} = readdlm(s_matrix_file, ',', Float64)
    
    sym_factor = sqrt.(ds ./ norm.(fv))

    for i in eachindex(Γ[:, 1])
        # @show size(Γ[:, i])
        Γ[i,i] -= dot(Γ[:, i], sym_factor) / sym_factor[i]
        #Γ[i,i] -= sum(Γ[:, i])
    end

    open(s_matrix_file, "w") do file
        writedlm(file, Γ, ",")
    end


    plot_matrix(Γ, arclengths)


end

main()