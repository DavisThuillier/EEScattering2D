using CurveFit
using CSV
using DataFrames
using Plots
using LaTeXStrings

function main()
    data_dir = joinpath(@__DIR__, "..", "data", "gamma", "200_30", "umklapp")
    filename = joinpath(data_dir, "resistivity.csv")
    data = CSV.read(filename, DataFrames.DataFrame, comment = "#")
    data.resistivity

    fit = curve_fit(LinearFit, data.T .^ 2, data.resistivity)
    polyfit = curve_fit(Polynomial, data.T, data.resistivity, 2)

    domain = LinRange(first(data.T), last(data.T), 100)

    plt = plot(data.T, data.resistivity, seriestype = :scatter, label = false)
    plot!(plt, domain, fit.(domain .^2), label = L"\rho = A + B T^2")
    plot!(plt, domain, polyfit.(domain), label = L"\rho = A + B T + C T^2")
    plot!(plt, xlabel = L"T/T_F", ylabel = L"\rho (\mathrm{\Omega \cdot m})")

    # plt2 = plot(data.T.^2, data.resistivity, seriestype = :scatter)
    # plot!(plt2, domain .^2, fit.(domain .^2))
    # plot!(plt2, xlabel = L"(T/T_F)^2", ylabel = L"\rho (\Omega \cdot \mathrm{m})")

    display(plt)
    savefig(plt, joinpath("plots","resistivity.pdf") )

    @show fit
    @show polyfit
    #println("ρ = ", fit[1], " + ", fit[2], " (T/T_F)^2")
end

main()