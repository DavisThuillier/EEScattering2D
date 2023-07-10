# [src/EEScattering2D.jl]
module EEScattering2D

    export helloworld, FermiSurfaceMesh

    import StaticArrays: SVector

    include("mesh.jl")
    using .FermiSurfaceMesh

    include("integration.jl")
    using .Integration

    function helloworld()
        println("I'm sorry, Dave. I'm afraid I can't do that.")
        return nothing
    end

end # module EEScattering2D
