# [src/EEScattering2D.jl]
module EEScattering2D

    export FermiSurfaceMesh #Integration
    export uniform_fermi_surface

    import StaticArrays: SVector
    using DelimitedFiles

    include("mesh.jl")
    import .FermiSurfaceMesh

    # include("integration.jl")
    # import .Integration

    function uniform_fermi_surface(thetas::Vector{Float64}, hamiltonian::Function, directory::String; write = false)
        fs = Vector{SVector{2, Float64}}(undef, length(thetas))
        fv = Vector{SVector{2, Float64}}(undef, length(thetas))

        FermiSurfaceMesh.fill_fermi_surface!(fs, thetas, hamiltonian)
        FermiSurfaceMesh.fill_fermi_velocity!(fv, fs, hamiltonian)

        write && begin
            open(joinpath(directory, "fermi_surface_$(length(fs)).csv"), "w") do file
                println(file, "kx,ky,dh/dx,dh/dy")
                writedlm(file, hcat(first.(fs), last.(fs), first.(fv), last.(fv)), ",")
            end
        end

        return fs, fv
    end
end # module EEScattering2D
