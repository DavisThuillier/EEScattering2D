# [src/EEScattering2D.jl]
module EEScattering2D

    export FermiSurfaceMesh, FermiSurfaceIntegration
    export uniform_fermi_surface

    import StaticArrays: SVector
    using DelimitedFiles

    include("mesh.jl")
    import .FermiSurfaceMesh

    include("integration.jl")
    import .FermiSurfaceIntegration

    # "Assuming injection angles between 0 and pi/4, construct the full collision matrix."
    # function generate_full_matrix(sub_mat::Matrix{Float64})
    #     num_columns = size(mat)[2]
    #     full_mat = Matrix{Float64}(undef, num_columns, num_columns)

    #     full_mat[]
    #     full_mat[2 + div(matrix_dim, 8): div(matrix_dim, 4) + 1, :] = circshift( reverse( n_matrix[1:div(matrix_dim, 8), :]), (0, 2 * div(matrix_dim, 8) + 1) )
    # end

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
