const band::String = "gamma"
const resolution::String = "75_25"
const matrix_dim::Int = 1200
const temperature::Float64 = 0.001581

# data_dir = joinpath(@__DIR__, "..", "..", "data", "$(band)_fermi_profile", "$(round(temperature, digits = 8))", resolution)
data_dir = joinpath(@__DIR__, "..", "..", "data","gamma_band", "$(temperature)", resolution)
ξ_matrix_file::String = joinpath(data_dir, "Γ_ξ_$(matrix_dim)_$(temperature).csv")
s_matrix_file::String = joinpath(data_dir, "Γ_s_$(matrix_dim)_$(temperature).csv")