# using LsqFit
using CurveFit
using CSV
using DataFrames
using Plots
using LaTeXStrings

fit_model(x::Vector{Float64}, p::Vector{Float64}) = x[2] - (p[1] + p[2] * x[1]^p[3] )

# fit_model(T::Vector{Float64}, p::Vector{Float64}) = p[1] * T .^ p[2] 

function main()
    filename = joinpath(@__DIR__, "resistivity.csv")
    data = CSV.read(filename, DataFrames.DataFrame, comment = "#")
    data.resistivity
    data.T *= 6326.35

    shift = 0
    logfit = curve_fit(LinearFit, log.(data.T), log.(data.resistivity .+ shift))
    square_fit = curve_fit(LinearFit, data.T.^2, data.resistivity)
    domain = collect(LinRange(first(data.T), last(data.T), 100))

    # border_ratio = 1.05

    plt = plot(data.T, data.resistivity, seriestype = :scatter, label = false)
    plot!(plt, domain, -shift .+ exp.(logfit.(log.(domain))))
    # plot!(plt, domain, square_fit.(domain.^2), label = "") #label = 
    plot!(plt, xlabel = L"T (K)", ylabel = L"\rho (\mathrm{n\Omega \cdot cm})")
    display(plt)

    # plt2 = plot(data.T.^2, data.resistivity, seriestype = :scatter, label = "")
    # #plot!(plt2, domain .^2, fit.(domain .^2))
    # plot!(plt2, xlabel = L"(T/T_F)^2", ylabel = L"\rho (\Omega \cdot \mathrm{m})", xlims = (0.0, border_ratio * last(domain)^2), ylims = (0.0, border_ratio * last(data.resistivity)))
    # display(plt2)

    
    # plt3 = plot(log.(data.T), log.(data.resistivity), seriestype = :scatter, legend = false)
    # plot!(plt3, xlabel = L"\ln(T/T_F)", ylabel = L"\ln(\rho/\rho_0)")
    # plot!(log.(domain), logfit.(log.(domain)) ,label = "")
    # display(plt3)

    #savefig(plt3, joinpath("plots","resistivity_log.pdf") )

    
    
    @show logfit
    @show square_fit
    # println("ρ = ", power_fit.param[1], " + ", power_fit[2], " (T/T_F)^2")
end

main()