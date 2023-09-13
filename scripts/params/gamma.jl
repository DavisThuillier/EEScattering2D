### Gamma Band ###

# In units of t0 = 119 meV
const mu::Float64 = 1.48 
const tx::Float64 = 1.0 
const ty::Float64 = 1.0 
const tp::Float64 = 0.392 
const ef::Float64 = 2.0*tx + 2.0*ty + 4*tp + mu
# T_f =  (ef * 0.65 / 8.617333e-5 * 0.119)# K

hamiltonian(k::SVector{2, Float64}) = - 2.0 * (tx/ef) * cos(2pi*k[1]) - 2.0 * (ty/ef) * cos(2pi * k[2]) - 4 * (tp/ef) * cos(2pi * k[1]) * cos(2pi * k[2]) - mu/ef

### Deformation Potentials ###
const alph::Float64 = 7.604
const alph_p::Float64 = 7.604

function dxx(k::SVector{2,Float64}, v::SVector{2,Float64})
    return 2 * alph * (tx/mu) * cos(2pi*k[1]) + 2 * alph_p * (tp/mu) * cos(2pi*k[1]) * cos(2pi*k[2]) - k[1] * v[1]
end

function dyy(k::SVector{2,Float64}, v::SVector{2,Float64})
    return 2 * alph * (ty/mu) * cos(2pi*k[2]) + 2 * alph_p * (tp/mu) * cos(2pi*k[1]) * cos(2pi*k[2]) - k[2] * v[2]
end

function dxy(k::SVector{2,Float64}, v::SVector{2,Float64})
    return - 2 * alph_p * (tp/mu) * sin(2pi*k[1]) * sin(2pi*k[2]) - k[1] * v[2]
end