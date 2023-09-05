using Plots

function main()
    n_t = [13, 17, 21, 25, 29]
    gamma = [0.0125, 0.0113, 0.0104, 0.0099, 0.0095]

    plt = plot(n_t, gamma)
    display(plt)
end

main()
