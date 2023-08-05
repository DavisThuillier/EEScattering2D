band::String = "free"
resolution::String = "150_25"
umklapp::Bool = false
umklapp_dir::String = "umklapp"
!umklapp && (umklapp_dir = "ohne_"*umklapp_dir)
matrix_dim::Int = 800
temperature::Float64 = 0.016
run::String = "2023.07.25(12)"

data_dir = joinpath(@__DIR__, "..", "..", "data", band, resolution, umklapp_dir, "$(round(temperature, digits = 4))", run)
mat_filename::String = joinpath(data_dir, "Γ_full_$(matrix_dim)_$(temperature).csv")