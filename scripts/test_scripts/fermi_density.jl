using Plots

fd(E::Float64) = 1 / (exp(E) + 1)

function main(N::Int)
    iseven(N) && (N += 1)

    α::Float64 = 0.001
    β::Float64 = sqrt(1 - 4*α)

    E_max::Float64 = 1

    energies = Vector{Float64}(undef, N)
    energies[1] = 0.0
    for i in 2:div(N,2)+1
        energies[i] = log(1 / (fd(energies[i - 1]) - β/N) - 1)
        energies[N - i + 2] = - energies[i] 
    end
    circshift!(energies, div(N,2))

    plot(energies, fd.(energies), seriestype = :scatter)
end

main(30)