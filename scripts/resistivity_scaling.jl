using CurveFit
using CSV
using DataFrames
using Plots
using LaTeXStrings

function main()
    data_dir = joinpath(@__DIR__, "..", "data", "gamma", "150_25", "umklapp")
    filename = joinpath(data_dir, "resistivity.csv")
    data = CSV.read(filename, DataFrames.DataFrame, comment = "#")
    data.resistivity

    power_fit = curve_fit(PowerFit, data.T, data.resistivity)
    fit = curve_fit(LinearFit, data.T .^ 2, data.resistivity)
    polyfit = curve_fit(Polynomial, data.T, data.resistivity, 2)

    domain = LinRange(first(data.T), last(data.T), 100)

    border_ratio = 1.05

    plt = plot(data.T, data.resistivity, seriestype = :scatter, label = false)
    plot!(plt, domain, power_fit.(domain))
    #plot!(plt, domain, fit.(domain .^2), label = L"\rho = A + B T^2")
    #plot!(plt, domain, polyfit.(domain), label = L"\rho = A + B T + C T^2")
    plot!(plt, xlabel = L"T/T_F", ylabel = L"\rho (\mathrm{\Omega \cdot m})")
    display(plt)

    plt2 = plot(data.T.^2, data.resistivity, seriestype = :scatter, label = "")
    #plot!(plt2, domain .^2, fit.(domain .^2))
    plot!(plt2, xlabel = L"(T/T_F)^2", ylabel = L"\rho (\Omega \cdot \mathrm{m})", xlims = (0.0, border_ratio * last(domain)^2), ylims = (0.0, border_ratio * last(data.resistivity)))
    display(plt2)

    plt3 = plot(log.(data.T), log.(data.resistivity), seriestype = :scatter)
    plot!(plt3, xlabel = L"\ln(T/T_F)", ylabel = L"\ln(\rho/\rho_0)")
    display(plt3)

    #savefig(plt3, joinpath("plots","resistivity_log.pdf") )

    logfit = curve_fit(LinearFit, log.(data.T), log.(data.resistivity))
    
    @show logfit
    @show power_fit
    #println("ρ = ", fit[1], " + ", fit[2], " (T/T_F)^2")
end

main()