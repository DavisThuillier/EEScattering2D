# [src/EEScattering2D.jl]
module EEScattering2D

    export FermiSurfaceMesh, FermiSurfaceIntegration
    export uniform_fermi_surface, get_dos, inner_product

    import StaticArrays: SVector
    import LinearAlgebra: norm
    using DelimitedFiles

    include("mesh.jl")
    import .FermiSurfaceMesh

    include("integration.jl")
    import .FermiSurfaceIntegration

    function inner_product(vec1, vec2, fs, hamiltonian::Function, T::Float64)
        ds = FermiSurfaceMesh.get_ds(fs)
        weights = FermiSurfaceIntegration.fd.(hamiltonian.(fs), T) .* (1 .- FermiSurfaceIntegration.fd.(hamiltonian.(fs), T))
        return vec1' * ((weights .* vec2) .* ds) 
    end

    function get_dos(fs::Vector{SVector{2,Float64}}, fv::Vector{SVector{2,Float64}})
        L::Int = length(fs)
        ds = FermiSurfaceMesh.get_ds(fs)
    
        density::Float64 = 0.0
        for i in eachindex(fs) 
            density += 0.5 * ds[i] * (1/norm(fv[i]) + 1/norm(fv[mod(i, L) + 1]))
        end
    
        return density
    end

    "Set diagonal elements to give null row sum asssuming uniform arclength grid."
    function enforce_particle_conservation!(mat::Matrix{Float64})
        for i in minimum(size(mat)) # Iterate over shorter matrix dimension
            mat[i,i] -= sum(matrix[i,:])
        end
    end

    "Set diagonal elements to give null row sum with nonuniform arclength grid."
    function enforce_particle_conservation!(mat::Matrix{Float64}, fs::Vector{SVector{Float64}})
        ds = FermiSurfaceMesh.get_ds(fs)
        for i in eachindex(fs) # Iterate over shorter matrix dimension
            mat[i,i] -= dot(matrix[i,:], ds) / ds[i]
        end
    end

end # module EEScattering2D
