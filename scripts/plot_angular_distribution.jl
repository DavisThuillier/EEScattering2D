using CSV
using DataFrames
using Plots
using Statistics
using LaTeXStrings
using FastTransforms
import CurveFit: linear_fit

const thermal_scale = cgrad([:blue, :red], [0.25,0.5,0.75]) # Color gradient for the plots

function nUFFT(angles, vals, m::Int64, central_angle)
    lambda_m = 0

    for i in 1:(length(angles) - 1)  # Naive NUFFT
            delta = angles[i + 1] - angles[i]
            lambda_m += vals[i] * cos( m * (angles[i] - central_angle)) * delta
    end

    return lambda_m / pi
end

function normalized_temp(temperatures, i)
    # For calculating color value
    return (log2(temperatures[i] / temperatures[1]) ) / (length(temperatures)  - 1)
end


function main()
    include("plot_params.jl")

    directory = joinpath("data", band, "$(num_angles)_$(p_num)", "theta_$(round(central_angle, digits = 4))")
    eigenvalues = Matrix{Complex}(undef, length(temperatures), length(modes))

    plt = plot(title = latexstring("\$\\sigma(\\theta), \\quad\\theta_0 = {$(round(central_angle, digits = 4))}\$"), legendtitle = L"T/T_F", fmt = :png)
    back_plt = plot(title = latexstring("\$\\sigma(\\theta) \times T / T_f\$"), legendtitle = L"T/T_F")

    for i in eachindex(temperatures)
        filename = joinpath(directory, "I_theta_$(round(temperatures[i], digits=4)).csv")
        data = CSV.read(filename, DataFrames.DataFrame)
        delete!(data, 1) # Remove the injection angle data point

        # Perform Fourier transform
        #println("T = $(temperatures[i]) T_F")
        gamma_0 = nUFFT(data.theta, data.I_theta, 0, central_angle) 
        for m in eachindex(modes)
            eigenvalues[i,m] = (nUFFT(data.theta, data.I_theta, modes[m], central_angle) - gamma_0) / temperatures[i]^2
            #println("lambda_$(modes[m]) = ", round(eigenvalues[i,m], digits = 4))
        end
        #println()

        ref_radius = median(abs.(eigenvalues[1, :])) /4 # Radius of the circle from which deviations indicate the angular scattering
        
        plot!(plt, data.theta, data.I_theta .+ ref_radius, labels = round(temperatures[i], digits = 4), color = get(thermal_scale, normalized_temp(temperatures,i)), legend = :outerbottomleft, proj = :polar, linewidth = 0.4, ylims = (0, 1.5 * ref_radius))

        plot!(back_plt, (data.theta .- (central_angle + pi)) / pi / temperatures[i], data.I_theta / temperatures[i], labels = round(temperatures[i], digits = 4), color = get(thermal_scale, normalized_temp(temperatures, i) ), xlims = (-2,2), xticks = (-3:3, 
        [L"-3\theta_T", L"-2\theta_T", L"-\theta_T", "0", L"\theta_T", L"2\theta_T", L"3\theta_T"]), ylims = (-40, 5))# Back_scattering region scaled in angle by temperature
    end

    plot!(plt, dpi = 500)
    savefig(plt, joinpath(directory, "sigma_$(num_angles)_$(p_num).png"))
    display(plt)
    #display(back_plt)

    ##################
    ### Linear Fit ###
    ##################

    leigenvalues = log.(-real.(eigenvalues[6:end, :]))
    ltemperatures = log.(temperatures[6:end])

    params3 = linear_fit(ltemperatures, leigenvalues[:,2])
    params5 = linear_fit(ltemperatures, leigenvalues[:,4])
    
    eigen_title::String = "Gamma"

    eigen_plot = plot(log.(temperatures), log.(-real.(eigenvalues)), labels = transpose(modes), legend =:outertopright, xlabel = latexstring("Temperature, \$\\ln(T/T_F)\$"), ylabel = latexstring("Eigenvalues, \$\\ln(\\lambda_m T^2_F/T^2)\$"), color = mode_colors, title = latexstring("$(eigen_title) \$N = $(num_angles) \\times $(p_num),\\theta_0 = $(round(central_angle, digits = 4))\$"), dpi = 500)

    fit_domain = collect(-4:0.5:-2.5)
    plot!(fit_domain, params3[1] .+ params3[2] .* fit_domain, color = :black, linestyle = :dash, label = latexstring("\$T^{$(round(params3[2] + 2, digits = 3))}\$"))
    plot!(fit_domain, params5[1] .+ params5[2] .* fit_domain, color = :blue, linestyle = :dash, label = latexstring("\$T^{$(round(params5[2] + 2, digits = 3))}\$"))

    savefig(eigen_plot, joinpath(directory, "eigens_$(num_angles)_$(p_num)"))
end

main()