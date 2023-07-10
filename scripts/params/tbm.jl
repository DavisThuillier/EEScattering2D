### Example Band ###

const mu::Float64 = 0.2
const tx::Float64 = 1.0
const ty::Float64 = 1.0

hamiltonian(k::SVector{2, Float64}) = - 2.0 * tx * cos(k[1]*pi) - 2.0 * ty * cos(k[2]*pi) - mu