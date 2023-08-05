### Free Electron Gas ###

const m::Float64 = 0.10335 # Effective mass approximates the gamma band of Sr2RuO4

hamiltonian(k::SVector{2, Float64}) = 0.5 * (k[1]^2 + k[2]^2) / m - 1.0