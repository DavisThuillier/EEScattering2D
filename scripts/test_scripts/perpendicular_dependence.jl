using EEScattering2D
using Plots
using StaticArrays
using LaTeXStrings

function main()
    n_t = 41
    n_s = 208

    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim)
    perimeter = FermiSurfaceMesh.get_perimeter(fs)
    fv = Vector{SVector{2, Float64}}(undef, length(fs))
    FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  

    ef::Float64 = 0.55 # Fermi energy in eV
    e0::Float64 = 55.26349406e6 # Vacuum permittivity in e^2 / eV / m
    c::Float64  = 12.68e-10 # z-axis dimension of unit cell in meters

    alpha = 1 / (ef * e0 * c) # Characteristic non-dimensionalized energy scale for interaction matrix element
    q_squared = alpha * get_dos(fs, fv) / (2pi)^2 # Thomas-Fermi screening wavevector squared 

    full_matrix = Matrix{Float64}(undef, n_t - 2, n_t - 2)

    s_1 = 0.0
    for s_2 in 1:100
        loci = [s_1, s_2 * perimeter / n_s]
        integration_mesh, mesh_dVs, mesh_variance, arclengths , loci_indices = FermiSurfaceMesh.discretize(fs, n_s, n_t, loci, hamiltonian, temperature, prec) 

        integral::SVector{2,Float64} = [0.0,0.0]
        for i in 2:(n_t-1)
            p1 = (integration_mesh[end,loci_indices[1]] - integration_mesh[begin, loci_indices[1]]) * i / n_t + integration_mesh[begin, loci_indices[1]]
            energies = hamiltonian.(energies)
            for j in 2:(n_t-1)
                k = (integration_mesh[end,loci_indices[2]] - integration_mesh[begin, loci_indices[2]]) * j / n_t + integration_mesh[begin, loci_indices[2]]
                integral = FermiSurfaceIntegration.collision_integral(p1, k, integration_mesh, energies, dVs, hamiltonian, variance, temperature, q_squared, umklapp = umklapp)
                full_matrix[i,j] = sum(integral)
            end
        end

        plt4 = heatmap(full_matrix, clims = (-1,4), cmap = :berlin, aspect_ratio = 1.0)
        plot!(plt4, xlabel = L"t_1", ylabel = L"t_2", title = latexstring("\$ s_1 = $(s_1), s_2 = $(s_2)\$"))
        display(plt4)
    end
    
    # FermiSurfaceIntegration.collision_integral((1, div(perp_num,2)), momenta, energies, dVs, hamiltonian, variance, temperature, q_squared, umklapp = umklapp)
        

end

include("../params/params.jl")
include(joinpath(@__DIR__, "..", "params", "$(band).jl"))

main()