# using LsqFit
using CurveFit
using CSV
using DataFrames
using Plots
using LaTeXStrings

fit_model(T::Vector{Float64}, p::Vector{Float64}) = p[1] .+ p[2] * T.^p[3] 

# fit_model(T::Vector{Float64}, p::Vector{Float64}) = p[1] * T .^ p[2] 

function main()
    filename = joinpath(@__DIR__, "resistivity.csv")
    data = CSV.read(filename, DataFrames.DataFrame, comment = "#")
    data.resistivity

    power_fit = curve_fit(PowerFit, data.T, data.resistivity)
    fit = curve_fit(LinearFit, data.T .^ 2, data.resistivity)
    # polyfit = curve_fit(Polynomial, data.T, data.resistivity, 2)

    # fit = curve_fit(fit_model, data.T, data.resistivity, [0.0, 1.0, 2.0])
    # @show fit.param
    domain = collect(LinRange(first(data.T), last(data.T), 100))

    # border_ratio = 1.05

    plt = plot(data.T, data.resistivity, seriestype = :scatter, label = false)
    # plot!(plt, domain, power_fit.(domain), label = latexstring("\$\\rho_{xx} = $(power_fit[1]) (T / T_F)^{$(power_fit[2])}\$"))
    plot!(plt, domain, power_fit.(domain), label = "") #label = latexstring("\$\\rho_{xx} = $(round(fit[1], digits = 4)) T^{$(round(fit[2], digits = 4))}\$"))
    #plot!(plt, domain, fit.(domain .^2), label = L"\rho = A + B T^2")
    #plot!(plt, domain, polyfit.(domain), label = L"\rho = A + B T + C T^2")
    plot!(plt, xlabel = L"T/T_F", ylabel = L"\rho (\mathrm{n\Omega \cdot cm})")
    display(plt)

    # plt2 = plot(data.T.^2, data.resistivity, seriestype = :scatter, label = "")
    # #plot!(plt2, domain .^2, fit.(domain .^2))
    # plot!(plt2, xlabel = L"(T/T_F)^2", ylabel = L"\rho (\Omega \cdot \mathrm{m})", xlims = (0.0, border_ratio * last(domain)^2), ylims = (0.0, border_ratio * last(data.resistivity)))
    # display(plt2)

    logfit = curve_fit(LinearFit, log.(data.T), log.(data.resistivity))
    plt3 = plot(log.(data.T), log.(data.resistivity), seriestype = :scatter, legend = false)
    plot!(plt3, xlabel = L"\ln(T/T_F)", ylabel = L"\ln(\rho/\rho_0)")
    plot!(log.(domain), logfit.(log.(domain)) ,label = "")
    display(plt3)

    #savefig(plt3, joinpath("plots","resistivity_log.pdf") )

    
    
    @show logfit
    @show power_fit
    println("ρ = ", fit[1], " + ", fit[2], " (T/T_F)^2")
end

main()