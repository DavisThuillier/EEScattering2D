using DelimitedFiles
using CSV
using DataFrames
using Plots
using CurveFit

function main()
    data_dir = joinpath(@__DIR__, "..", "data")
    matrix_dim = 1200
    res = "192_41"

    Tf = 6326.35178
    temperatures = collect(8.0:0.25:12) / Tf
    min_widths = zeros(Float64, length(temperatures))
    for i in eachindex(temperatures)
        temp_string = string(round(temperatures[i], digits = 6))
        filename = joinpath(data_dir, temp_string, res, "Γ_full_$(matrix_dim)_"*temp_string*".csv")
        
        if isfile(filename)
            full_matrix = readdlm(filename,',',Float64)
            width = Vector{Float64}(undef, matrix_dim)
            for i in 1:size(full_matrix)[1]
                max, index = findmax(full_matrix[i, :])
                j = index
                while full_matrix[j] < max / 2.0
                    # j = mod(j, matrix_dim) + 1 
                    j = mod(j - 2, matrix_dim) + 1                   
                end
                width[i] = (j - 1 + (max / 2.0 - full_matrix[j - 1]) /(full_matrix[j] - full_matrix[j - 1])) 
            end
            # plt = plot(width, ylims = (0.0,10.0), title = temperatures[i])
            # display(plt)
            min_widths[i] = minimum(width)
        end
    end

    outfile = joinpath(@__DIR__, "peak_widths.csv")
    open(outfile, "w") do file
        println(file,"T,hwhm")
        writedlm(file,hcat(temperatures,min_widths), ",")
    end

    plt = plot(temperatures, min_widths, seriestype = :scatter)
    display(plt)
end

function fit_widths()
    infile = joinpath(@__DIR__, "peak_widths.csv")
    data = CSV.read(infile, DataFrames.DataFrame)
    
    fit = curve_fit(LinearFit, data.T, data.hwhm)
    @show fit

    plt = plot(data.T, fit.(data.T))
    plot!(plt, data.T, data.hwhm, seriestype = :scatter)
end

# main()
fit_widths()