using EEScattering2D
using Plots
using StaticArrays
using LaTeXStrings

function main()
    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim)
    fv = Vector{SVector{2, Float64}}(undef, length(fs))
    FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)  

    ef::Float64 = 0.55 # Fermi energy in eV
    e0::Float64 = 55.26349406e6 # Vacuum permittivity in e^2 / eV / m
    c::Float64  = 12.68e-10 # z-axis dimension of unit cell in meters

    alpha = 1 / (ef * e0 * c) # Characteristic non-dimensionalized energy scale for interaction matrix element
    prefactor = alpha^2 
    q_squared = alpha * get_dos(fs, fv) / (2pi)^2 # Thomas-Fermi screening wavevector squared 

    momenta, dVs, variance, arclengths = FermiSurfaceMesh.discretize(fs, num_bins, perp_num, 1, hamiltonian, temperature, prec)

    plt = plot(first.(momenta), last.(momenta), seriestype= :scatter, markersize= 0.05, legend = false, markershape = :diamond, aspect_ratio = :equal, title = "T = $(round(temperature, digits = 4))", xlims = (-0.75,0.75), ylims = (-0.75,0.75))
    plot!(plt, first.(fs), last.(fs), color = :blue)
    display(plt)

    energies = hamiltonian.(momenta)

    full_matrix = Matrix{Float64}(undef, perp_num, perp_num)
    df_matrix   = Matrix{Float64}(undef, size(momenta))

    t_1 = 2
    
    # for i in CartesianIndices(momenta)
    #     k_index = (i[1], i[2])
    #     integral = FermiSurfaceIntegration.collision_integral((t_1, s_1), k_index, momenta, energies, dVs, hamiltonian, variance, temperature, q_squared, umklapp = umklapp)
    #     full_matrix[i[1], i[2]] = sum(integral)
    #     df_matrix[i[1], i[2]] = 1 / FermiSurfaceIntegration.fd_normalization(hamiltonian(momenta[i[1], i[2]]), temperature)
    # end

    # @show df_matrix[1, 1] / df_matrix[div(perp_num, 2), 1]

    # plt = scatter(first.(momenta), last.(momenta), zcolor = full_matrix, markersize = 1, markerstrokewidth = 0.0, markershape = :diamond)

    # plt = heatmap(full_matrix, cmap = :berlin)
    # #plt = surface(first.(momenta), last.(momenta), full_matrix)
    # plot!(plt, xlabel = L"s_k", ylabel = L"t_k")
    # display(plt)

    # plt2 = heatmap(df_matrix, cmap = :berlin)
    # display(plt2)

    s_1 = 1
    for s_2 in 1:num_bins

        integral::SVector{2,Float64} = [0.0,0.0]
        for i in eachindex(momenta[:, s_1])
            for j in eachindex(momenta[:, s_2])
                integral = FermiSurfaceIntegration.collision_integral((i, s_1), (j, s_2), momenta, energies, dVs, hamiltonian, variance, temperature, q_squared, umklapp = umklapp)
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