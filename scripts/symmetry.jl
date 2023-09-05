using EEScattering2D

import StaticArrays: SVector
using LaTeXStrings
using Plots
using ProgressBars
import LinearAlgebra: norm
using DelimitedFiles

function main()
    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim)
    arclengths = FermiSurfaceMesh.get_arclengths(fs)
    fv = Vector{SVector{2, Float64}}(undef, length(fs))
    FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)

    perimeter = FermiSurfaceMesh.get_perimeter(fs)

    momenta, dVs , var , arclengths, _ = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, 1, hamiltonian, temperature, prec)

    ef::Float64 = 0.55 # Fermi energy in eV
    e0::Float64 = 55.26349406e6 # Vacuum permittivity in e^2 / eV / m
    c::Float64  = 12.68e-10 # z-axis dimension of unit cell in meters

    alpha = 1 / (ef * e0 * c) # Characteristic non-dimensionalized energy scale for interaction matrix element
    q_squared = alpha * get_dos(fs, fv) / (2pi)^2 # Thomas-Fermi screening wavevector squared 

    n_s = size(momenta)[2]
    n_t = perp_num

    asymmetry = Matrix{Float64}(undef, n_s, n_s)

    for i in ProgressBar(1:n_s)
        s1 = arclengths[i]
        for j in 1:n_s
            s2 = arclengths[j]
            integration_mesh, mesh_dVs, mesh_variance, _ , loci_indices = FermiSurfaceMesh.discretize(fs, n_s, n_t, [s1, s2], hamiltonian, temperature, prec)

            k = integration_mesh[div(n_t, 2) + 1, loci_indices[2]]

            C1 = [0.0,0.0]
            C2 = [0.0,0.0]

            for m in 2:(n_t - 1)
                p1 = integration_mesh[m, loci_indices[1]]
                dt1 = norm(integration_mesh[m + 1, loci_indices[1]] - integration_mesh[m - 1, loci_indices[1]]) / 2
                for n in 2:(n_t - 1)
                    k = integration_mesh[n, loci_indices[2]]
                    dtk = norm(integration_mesh[n + 1, loci_indices[2]] - integration_mesh[n - 1, loci_indices[2]]) / 2
                    C1 += FermiSurfaceIntegration.collision_integral(p1, k, integration_mesh, hamiltonian.(integration_mesh), mesh_dVs, hamiltonian, mesh_variance, temperature, q_squared) * dt1 * dtk
                    C2 += FermiSurfaceIntegration.collision_integral(k, p1, integration_mesh, hamiltonian.(integration_mesh), mesh_dVs, hamiltonian, mesh_variance, temperature, q_squared) * dt1 * dtk
                end
            end
            asymmetry[i,j] = abs(sum(C1 - C2)) / abs(sum(C1 + C2))
            println("Asymmetry at s1 = $(round(s1,digits = 4)), s2 = $(round(s2,digits=4)):  ", round(asymmetry[i,j],digits=6))
        end
        plt = plot(arclengths / perimeter, asymmetry[i,:], seriestype = :scatter)
        plot!(plt, title = latexstring("Asymmetry, \$s_1 = $(s1)\$"))
        display(plt)
    end

    data_dir = joinpath(@__DIR__, "..", "data", "$(band)_fermi_profile", "$(temperature)", "$(num_bins)_$(perp_num)")
    !isdir(data_dir) && mkpath(data_dir)
    outfile = "asymmetry.csv"
    open(joinpath(data_dir, outfile), "w") do file
        writedlm(file, asymmetry, ",")
    end
end

function plot_asymmetry()
    data_dir = joinpath(@__DIR__, "..", "data", "$(band)_fermi_profile", "$(temperature)", "$(num_bins)_$(perp_num)")
    infile = joinpath(data_dir,"asymmetry.csv")

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim)
    perimeter = FermiSurfaceMesh.get_perimeter(fs)
    _, _ , _ , arclengths, _ = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, 1, hamiltonian, temperature, prec)
    
    asymmetry_matrix::Matrix{Float64} = readdlm(infile, ',', Float64)
    plt = heatmap(arclengths / perimeter, arclengths / perimeter, asymmetry_matrix, cmap = :berlin, aspect_ratio = 1.0)
    plot!(plt, xlabel = L"s_1", ylabel = L"s_2", title = "Asymmetry")
    display(plt)
end

include("params/params.jl")
include("params/$(band).jl")


# main()
plot_asymmetry()