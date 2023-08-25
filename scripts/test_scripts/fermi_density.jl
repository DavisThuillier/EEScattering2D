using Plots

fd(E::Float64) = 1 / (exp(E) + 1)

function main(Ns::Vector{Int})
    for N in Ns
        iseven(N) && (N += 1)

        α::Float64 = 0.01
        β::Float64 = sqrt(1 - 4*α)

        energies = Vector{Float64}(undef, N)
        energies[1] = 0.0
        for i in 2:div(N,2)+1
            energies[i] = log(1 / (fd(energies[i - 1]) - β/(N-1)) - 1)
            energies[N - i + 2] = - energies[i] 
        end
        circshift!(energies, div(N,2))

        # integral = 0.0
        # for i in 2:N-1
        #     integral += fd(energies[i]) * (1 - fd(energies[i])) * (energies[i + 1] - energies[i - 1])/2
        # end
        # @show integral
        e_bound = (energies[2] + energies[1]) / 2
        @show fd(e_bound) - fd(-e_bound)

        int_domain = LinRange(-e_bound, e_bound, 100)
        full_domain = LinRange(first(energies), last(energies), 100)
        plt = plot(energies, fd.(energies), seriestype = :scatter, title = "N = $(N)")
        plot!(plt, full_domain, fd.(full_domain), linestyle = :dash)
        plot!(plt, int_domain, fd.(int_domain), linewidth = 2.0)
        display(plt)
    end
end

main([9,13,17,21,25,29])