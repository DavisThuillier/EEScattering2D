const band::String = "gamma"
const resolution::String = "120_21"
const umklapp::Bool = true
const umklapp_dir::String = "umklapp"
!umklapp && (umklapp_dir = "ohne_"*umklapp_dir)
const matrix_dim::Int = 1200
const temperature::Float64 = 0.001581

# data_dir = joinpath(@__DIR__, "..", "..", "data", "$(band)_fermi_profile", "$(round(temperature, digits = 8))", resolution)
data_dir = joinpath(@__DIR__, "..", "..", "data","gamma_band", "$(temperature)", resolution)
mat_filename::String = joinpath(data_dir, "Γ_full_$(matrix_dim)_$(temperature).csv")