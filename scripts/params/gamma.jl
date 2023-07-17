### Gamma Band ###

# In units of t0 = 119 meV
const mu::Float64 = 1.48
const tx::Float64 = 1.0
const ty::Float64 = 1.0
const tp::Float64 = 0.41

hamiltonian(k::SVector{2, Float64}) = - 2.0 * (tx/mu) * cos(k[1]*pi) - 2.0 * (ty/mu) * cos(k[2]*pi) - 4 * (tp/mu) * cos(k[1]*pi) * cos(k[2]*pi)  - 1.0
