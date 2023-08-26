using Plots
using Roots
import SpecialFunctions: erf

function secant_root(f::Function, x0::Float64, x1::Float64, iterations::Int)
    x2::Float64 = 0.0
    for i in 1:iterations
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        x0, x1 = x1, x2
        abs(x0 - x1) < 1e-8 && break
    end 
    return x2
end

function main(loci::Vector{Float64}, perimeter::Float64)
    loci = mod.(vcat(loci, loci .+ perimeter / 2), perimeter)
    sort!(loci)
    
    walls = Vector{Float64}(undef, length(loci))
    for i in eachindex(loci)
        i == length(loci) && continue
        walls[i] = mod((loci[mod(i, length(loci)) + 1] + loci[i]) / 2, perimeter)
    end
    walls[end] = mod(((loci[1] + perimeter) + loci[end]) / 2, perimeter)

    regions = Vector{Tuple{Float64, Float64}}(undef, 2*length(loci))
    for i in eachindex(regions)
        if isodd(i)
            regions[i] = (loci[div(i, 2) + 1], walls[div(i, 2) + 1])
        else
            regions[i] = (walls[div(i,2)], loci[mod(div(i, 2), length(loci)) + 1])
        end
    end
    regions[end] = (regions[end][1], regions[end][2] + perimeter)

    sigma = 0.02
    A = 0.2
    N = 100
    mesh = Vector{Float64}(undef, 0)
    for i in eachindex(regions)
        locus = (i != length(regions)) ? loci[div(i, 2) + 1] : loci[1] + perimeter
        s = locus
        j = 1
        isapprox(regions[i][2], regions[i][1]) && continue  

        push!(mesh, regions[i][1])
        if isodd(i)
            while true
                s = secant_root(x -> ((x - locus) + A * erf((x - locus) / sigma)) - j/100, regions[i][1], (regions[i][2] + regions[i][1]) / 2, N)
                s > regions[i][2] && break
                push!(mesh, s)
                j += 1
            end
        else
            while true
                s = secant_root(x -> ((x - locus) + A * erf((x - locus) / sigma)) + j/100, regions[i][1], regions[i][2], N)
                regions[i][1] > s && break 
                push!(mesh, s)
                j += 1
            end
        end
    end
    sort!(mesh)

    domain = LinRange(0.0, 2*pi, 100)
    plt = plot(domain, map(x -> 1, domain), proj = :polar)
    plot!(plt, map(x -> (x / perimeter * 2pi, 1), loci), seriestype = :scatter)
    for i in 1:(length(regions) - 1)
        plot!(plt, LinRange(regions[i][1], regions[i][2] , 100) ./ (perimeter / 2pi), ones(100) )
    end
    plot!(plt, map(x -> (x / perimeter * 2pi, 1), mesh), seriestype = scatter, markershape = :cross, markersize = 0.5, color = :black)
    display(plt)
    # plot!(plt, map(x -> (x / perimeter * 2pi, 1), walls), seriestype = :scatter)
    
end

for i in range(0.0, 1.0, 21)
    for j in range(0.5, 1.5, 21)
        main([0.0, i, j],1.0)
    end
end