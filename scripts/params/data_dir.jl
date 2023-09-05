const band::String = "gamma"
const resolution::String = "128_33"
const umklapp::Bool = true
const umklapp_dir::String = "umklapp"
!umklapp && (umklapp_dir = "ohne_"*umklapp_dir)
const matrix_dim::Int = 1000
const temperature::Float64 = 0.002
const run::String = "2023.08.24"

# data_dir = joinpath(@__DIR__, "..", "..", "data", "$(band)_interp_test", resolution, umklapp_dir, "$(round(temperature, digits = 4))", run)
# data_dir = joinpath(@__DIR__, "..", "..", "data", "$(band)", "$(round(temperature, digits = 4))", resolution)
data_dir = joinpath(@__DIR__, "..", "..", "data", "$(band)_fermi_profile", "$(round(temperature, digits = 4))", resolution)
mat_filename::String = joinpath(data_dir, "Γ_full_$(matrix_dim)_$(temperature).csv")