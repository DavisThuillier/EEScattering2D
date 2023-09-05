using EEScattering2D

import StaticArrays: SVector
import LinearAlgebra: norm
using Plots

function main(s1::Float64, s2::Float64)
    fs = FermiSurfaceMesh.generate_fermi_surface(hamiltonian, row_dim)
    arclengths = FermiSurfaceMesh.get_arclengths(fs)
    fv = Vector{SVector{2, Float64}}(undef, length(fs))
    FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)

    p1 = FermiSurfaceMesh.get_momentum(fs, arclengths, s1)
    k  = FermiSurfaceMesh.get_momentum(fs, arclengths, s2)

    perimeter = FermiSurfaceMesh.get_perimeter(fs)
    loci = [s1, s2] * perimeter

    ef::Float64 = 0.55 # Fermi energy in eV
    e0::Float64 = 55.26349406e6 # Vacuum permittivity in e^2 / eV / m
    c::Float64  = 12.68e-10 # z-axis dimension of unit cell in meters

    alpha = 1 / (ef * e0 * c) # Characteristic non-dimensionalized energy scale for interaction matrix element
    q_squared = alpha * get_dos(fs, fv) / (2pi)^2 # Thomas-Fermi screening wavevector squared 

    for n_s in 208:4:244
        println("n_s: ", n_s)
        for n_t in 17:4:41
            integration_mesh, mesh_dVs, mesh_variance, arclengths , loci_indices = FermiSurfaceMesh.discretize(fs, n_s, n_t, loci, hamiltonian, temperature, prec) 
            
            C = [0.0, 0.0]
            for i in 2:(n_t - 1)
                p1 = integration_mesh[i, loci_indices[1]]
                dp1 = norm(integration_mesh[i + 1, loci_indices[1]] - integration_mesh[i - 1, loci_indices[1]]) / 2
                for j in 2:(n_t - 1)
                    p1 = integration_mesh[j, loci_indices[2]]
                    dk = norm(integration_mesh[j + 1, loci_indices[2]] - integration_mesh[j - 1, loci_indices[2]]) / 2
                    C += FermiSurfaceIntegration.collision_integral(p1, k, integration_mesh, hamiltonian.(integration_mesh), mesh_dVs, hamiltonian, mesh_variance, temperature, q_squared) * dk * dp1
                end
            end
            @show C 
        end
    end

end 

include("params/params.jl")
include("params/$(band).jl")

main(0.25, 0.8)