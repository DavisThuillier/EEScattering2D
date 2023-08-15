using LinearAlgebra
using DelimitedFiles
using Plots
using CurveFit

function main()
    
    matrix_dim = 1000

    data_dir = joinpath(@__DIR__, "..", "data", "gamma", "150_25", "umklapp") 
    date1 = "2023.08.07"
    date2 = "2023.08.08"
    temps = [0.001, 0.002, 0.003, 0.004]
    runs = joinpath.(string.(temps), [date1, date2, date2, date2])

    # data_dir = joinpath(@__DIR__, "..", "data", "gamma", "200_30", "umklapp") 
    # date = "2023.08.05"
    # temps = [0.001, 0.002, 0.004, 0.008, 0.016]
    # runs = joinpath.(string.(temps), fill(date, length(temps)))

    eigennorms = Vector{Float64}(undef, length(runs))
    for i in eachindex(runs)
        mat_filename = joinpath(data_dir, runs[i], "Γ_full_$(matrix_dim)_$(temps[i]).csv")
        full_matrix::Matrix{Float64} = readdlm(mat_filename, ',', Float64)

        lambdas = reverse(eigvals(full_matrix))
        eigennorms[i] = norm(lambdas)
    end
    @show eigennorms

    fit = curve_fit(LinearFit, log.(temps), log.(eigennorms))
    domain = LinRange(first(temps), last(temps), 100)
    @show fit
    plt = plot(log.(temps), log.(eigennorms), seriestype = :scatter)
    plot!(plt, log.(domain), fit.(log.(domain)))
end

main()

