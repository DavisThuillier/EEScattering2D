band::String = "gamma"
resolution::String = "200_30"
umklapp::Bool = true
umklapp_dir::String = "umklapp"
!umklapp && (umklapp_dir = "ohne_"*umklapp_dir)
matrix_dim::Int = 1000
temperature::Float64 = 0.001
run::String = "2023.08.05"

data_dir = joinpath(@__DIR__, "..", "..", "data", band, resolution, umklapp_dir, "$(round(temperature, digits = 4))", run)
mat_filename::String = joinpath(data_dir, "Γ_full_$(matrix_dim)_$(temperature).csv")