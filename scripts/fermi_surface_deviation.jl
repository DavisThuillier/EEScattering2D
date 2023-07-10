using Plots
using LaTeXStrings
using DataFrames
using CSV
using Interpolations

function normalized_temp(temperatures, i)
    # For calculating color value
    return 2 * (log2(temperatures[i] / temperatures[1]) ) / (length(temperatures)  - 1)
end

const thermal_scale = cgrad([:blue, :red], [0.25,0.5,0.75]) # Color gradient for the plots
const temperatures = 0.01 * (2.0 .^ (0:0.5:1.5))
const scale::Float64 = 0.0015

function main()
    include("params.jl")
    data_dir = joinpath(@__DIR__,"data", band)
    data_dir = joinpath(data_dir, "$(num_angles)_$(p_num)")

    fs_filename = joinpath(data_dir, "fermi_surface_1000.csv")
    fs = CSV.read(fs_filename, DataFrames.DataFrame)

    angles = mod.(atan.(fs.ky, fs.kx), 2*pi)

    plt_n = plot(fs.kx, fs.ky, aspect_ratio = 1.0, label = "Fermi Surface", color = :black, legend = :topleft, xlabel = L"\pi/a_x", ylabel = L"\pi/a_y", dpi = 400, title = "Normal Scattering")
    plt_u = plot(fs.kx, fs.ky, aspect_ratio = 1.0, label = "Fermi Surface", color = :black, legend = :topleft, xlabel = L"\pi/a_x", ylabel = L"\pi/a_y", dpi = 400, title = "Umklapp Scattering")
    plt = plot(fs.kx, fs.ky, aspect_ratio = 1.0, label = "Fermi Surface", color = :black, legend = :topleft, xlabel = L"\pi/a_x", ylabel = L"\pi/a_y", dpi = 400, title = "All Events")

    plts = [plt_n, plt_u, plt]

    angle_dir::String = "theta_0.5262"
    for i in eachindex(temperatures)
        filename = joinpath(data_dir, angle_dir, "T_$(round(temperatures[i], digits = 4)).csv")
        sigma = CSV.read(filename, DataFrames.DataFrame)
        sort!(sigma)
        @show sigma 
        # append!(sigma, first(sigma, 1))
        # sigma[end,1] = sigma[end,1] + 2*pi
        
        itp_n = extrapolate(interpolate(sigma.theta, sigma.I_n, FritschCarlsonMonotonicInterpolation()), Flat())
        itp_u = extrapolate(interpolate(sigma.theta, sigma.I_u, FritschCarlsonMonotonicInterpolation()), Flat())
        itp = extrapolate(interpolate(sigma.theta, sigma.I_n .+ sigma.I_u, FritschCarlsonMonotonicInterpolation()), Flat())

        itps = [itp_n, itp_u, itp]
        for j in eachindex(itps)
            devx = fs.kx .+ (scale * itps[j].(angles) / temperatures[i]) .* cos.(angles)
            devy = fs.ky .+ (scale * itps[j].(angles) / temperatures[i]) .* sin.(angles)

            plot!(plts[j], devx, devy, color = get(thermal_scale, normalized_temp(temperatures, i)), label = latexstring("\$ T = $(round(temperatures[i], digits = 4)) \$"))
        end
    end

    display(plt_n)
    display(plt_u)
    display(plt)
    #savefig(plt, joinpath(data_dir, angle_dir, "fs_deviation.png"))
end 

main()
