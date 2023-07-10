### Gamma Band ###

const mu::Float64 = 1.48
const tx::Float64 = 1.0
const ty::Float64 = 1.0
const tp::Float64 = 0.41

hamiltonian(k::SVector{2, Float64}) = - 2.0 * tx * cos(k[1]*pi) - 2.0 * ty * cos(k[2]*pi) - 4 * tp * cos(k[1]*pi) * cos(k[2]*pi)  - mu
