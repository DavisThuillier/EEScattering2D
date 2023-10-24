const band_folder::String = "com_convergence_test"
const resolution::String = "192x41_145"
const matrix_dim::Int = 1600
const temperature::Float64 = 0.001897

# data_dir = joinpath(@__DIR__, "..", "..", "data", "$(band)_fermi_profile", "$(round(temperature, digits = 8))", resolution)
data_dir = joinpath(@__DIR__, "..", "..", "data",band_folder, "$(temperature)", resolution)
ξ_matrix_file::String = joinpath(data_dir, "Γ_ξ_$(matrix_dim)_$(temperature).csv")
s_matrix_file::String = joinpath(data_dir, "Γ_s_$(matrix_dim)_$(temperature).csv")