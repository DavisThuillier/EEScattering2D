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

end # module EEScattering2D
