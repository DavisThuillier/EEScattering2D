band::String = "gamma"
resolution::String = "100_9"
umklapp::Bool = true
umklapp_dir::String = "umklapp"
!umklapp && (umklapp_dir = "ohne_"*umklapp_dir)
matrix_dim::Int = 1000
const temperature::Float64 = 0.002
run::String = "2023.08.24"

data_dir = joinpath(@__DIR__, "..", "..", "data", "$(band)_interp_test", resolution, umklapp_dir, "$(round(temperature, digits = 4))", run)
mat_filename::String = joinpath(data_dir, "Γ_full_$(matrix_dim)_$(temperature).csv")