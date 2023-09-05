using VoronoiCells
using Plots
using Random
using GeometryBasics


function main()
    rng = Random.MersenneTwister()
    rect = Rectangle(Point2(-0.5,-0.5), Point2(0.5,0.5))

    n = 900
    points = [Point2(rand(rng), rand(rng)) .- 0.5 for _ in 1:n]

    @time tess = voronoicells(points, rect)
    @time areas = voronoiarea(tess)
    # scatter(points, markersize = 6, label = "generators")
    # annotate!([(points[n][1] + 0.02, points[n][2] + 0.03, Plots.text(n)) for n in 1:10])
    plot(tess, legend = :topleft) 

end

main()