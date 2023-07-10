using DelimitedFiles
using Plots
using LinearAlgebra
using DataFrames
using CSV
using LaTeXStrings

function mirror_symmetrize(basis::Vector{Vector{ComplexF64}}, mirror_index::Int)
    mirror_mat = Matrix{ComplexF64}(undef, 2, 2)
    for i in 1:2
        for j in 1:2
            mirror_mat[i,j] = dot(basis[i], circshift(reverse(basis[j]), 2 * mirror_index - 1)) 
        end
    end
    weights = eigvecs(mirror_mat)

    w1 = weights[1,1] * basis[1] + weights[2,1] * basis[2]
    w2 = weights[1,2] * basis[1] + weights[2,2] * basis[2]

    theta1 = mean(mod.(angle.(w1), pi))
    theta2 = mean(mod.(angle.(w2), pi))

    w1 = real.(exp(- im * theta1) * w1)
    w2 = real.(exp(- im * theta2) * w2)

    return w1 / norm(w1), w2 / norm(w2) 
end

function fft(v::Vector{Float64}, m::Int)
    int::ComplexF64 = 0.0
    step::Float64 = 2pi / length(v)

    for k in eachindex(v)
        int += exp( - im * (m * step * (k - 1)) ) * step * v[k]
    end
    return int
end

function main()
    band = "gamma"
    resolution = "200_30"
    matrix_dim = 800
    temperature = 0.04

    dtheta = 2*pi / matrix_dim
    thetas = collect(range(0.0, 2 * pi, step = dtheta))
    
    data_dir = joinpath(@__DIR__, "data", band, resolution)
    fs = CSV.read(joinpath(data_dir, "fermi_surface_$(matrix_dim).csv"), DataFrames.DataFrame)
    fs_norms = sqrt.(fs.kx .^ 2 + fs.ky .^2)

    n_stem = joinpath(data_dir, "Γn_$(matrix_dim)_$(temperature)")
    u_stem = joinpath(data_dir, "Γu_$(matrix_dim)_$(temperature)")

    n_files = String[]
    u_files = String[]
    for file in readdir(data_dir; join=true)
        startswith(file, n_stem) && push!(n_files, file)
        startswith(file, u_stem) && push!(u_files, file)
    end
    
    if (length(n_files) != length(u_files) ) 
        println("Error: Missing data")
        exit()
    end

    full_matrix = Matrix{Float64}(undef, matrix_dim, matrix_dim)
    n_matrix = Matrix{Float64}(undef, matrix_dim, matrix_dim)
    u_matrix = Matrix{Float64}(undef, matrix_dim, matrix_dim)
    for i in eachindex(n_files)
        start_index = parse(Int, split(basename(n_files[i]), "_")[4])
        end_index   = parse(Int, split(basename(n_files[i])[begin:(end-4)], "_")[6])
        
        n_matrix[start_index:end_index, :] = readdlm(n_files[i], ',', Float64)
        u_matrix[start_index:end_index, :] = readdlm(u_files[i], ',', Float64)
    end

    # Enforce particle conservation by ensuring row sum is null
    for i in 1:(div(matrix_dim, 8) + 1)
        n_matrix[i,i] -= sum(n_matrix[i,:])
        u_matrix[i,i] -= sum(u_matrix[i,:])
    end

    ### Generate Full Matrix ###
    n_matrix[2 + div(matrix_dim, 8): div(matrix_dim, 4) + 1, :] = circshift( reverse( n_matrix[1:div(matrix_dim, 8), :]), (0, 2 * div(matrix_dim, 8) + 1) )
    u_matrix[2 + div(matrix_dim, 8): div(matrix_dim, 4) + 1, :] = circshift( reverse( u_matrix[1:div(matrix_dim, 8), :]), (0, 2 * div(matrix_dim, 8) + 1) )

    n_matrix[1 + div(matrix_dim, 4) : div(matrix_dim, 2), :] = circshift( n_matrix[1 : div(matrix_dim, 4), :], (0, div(matrix_dim, 4)))
    u_matrix[1 + div(matrix_dim, 4) : div(matrix_dim, 2), :] = circshift( u_matrix[1 : div(matrix_dim, 4), :], (0, div(matrix_dim, 4)))

    n_matrix[1 + div(matrix_dim, 2) : matrix_dim, :] = circshift( n_matrix[1 : div(matrix_dim, 2), :], (0, div(matrix_dim, 2)))
    u_matrix[1 + div(matrix_dim, 2) : matrix_dim, :] = circshift( u_matrix[1 : div(matrix_dim, 2), :], (0, div(matrix_dim, 2)))


    u_matrix *= (2pi/ matrix_dim)
    n_matrix *= (2pi/ matrix_dim)
    full_matrix = ( n_matrix + u_matrix) # Prefactor corresponds to the differential in the angular integral when taking the product of full_matrix and a vector
    
    lambdas = reverse(eigvals(full_matrix))
    eigenvecs = reverse(eigvecs(full_matrix), dims = 2) # Order eigenvectors from smallest to largest
    display(plot(real.(lambdas[1:20]), color = :green, seriestype = :scatter))


    #@show eigenvecs[:, end]

    # plt = plot(title = latexstring("Eigenmodes of the \$\\gamma\$ Band 1600") )
    # count = 0
    scale = 4.0

    println("### Mode Analysis ###")
    for i in 2:2:20
        w1, w2 = mirror_symmetrize([eigenvecs[:, i], eigenvecs[:, i + 1]], div(matrix_dim, 4) + 1)

        # @show lambdas[2 * i + 1]
        # @show (full_matrix * w1 ./ w1)[1]

        # @show dot(w1, full_matrix * w1)
        # @show lambdas[2 * i + 1]

        normal_contribution = dot(w1, n_matrix * w1) / dot(w1, full_matrix * w1)
        umklapp_contribution = dot(w1, u_matrix * w1) / dot(w1, full_matrix * w1)
    
        println("n = ", i/2)
        @show normal_contribution
        @show umklapp_contribution

        
        # println("Fourier Components of w1:")
        # for j in 0:20
        #     res = fft(w1, j)
        #     println("m = ", j, ": ", round(res, digits = 6))
        # end

        # println("Fourier Components of w2:")
        # for j in 0:10
        #     res = fft(w2, j)
        #     println("m = ", j, ": ", round(res, digits = 6))
        # end
        # println()


        plt = plot(thetas, fs_norms .+ scale * w2, title = latexstring("\$ \\lambda = $(i/2) \$"))
        plot!(plt, thetas, fs_norms, color = :black)
        # plot!(plt, thetas,fs_norms .- sin.((2 * i - 1) * thetas))
        display(plt)

        # plt2 = plot(thetas, fs_norms .+ scale * real.(w1), title = "m = $(i)", proj = :polar)
        # plot!(plt2, thetas, fs_norms, color = :black)
        # # plot!(plt2, thetas, fs_norms .+ cos.((2 * i - 1) * thetas))
        # display(plt2)
    end

    # plt2 = plot(thetas, fs_norms .+ scale * real.(eigenvecs[:,1598]), proj = :polar)
    # plot!(plt2, thetas, fs_norms .+ scale * imag.(eigenvecs[:,1598]))
    # plot!(plt2, thetas, fs_norms)
    # display(plt2)

    # plt3 = plot(thetas, fs_norms .+ 0.25 * cos.(thetas), proj = :polar)
    # display(plt3)
    # for i in matrix_dim:-1:(matrix_dim - 5)
        # println("Eigenvalue ", 1600 - i, ": ", eigenvalues[i])
        

        # plot!(plt, thetas, fs_norms + scale * eigenvecs[:, i], proj = :polar, label = "$(count)")
            # count += 1
        # if abs(eigenvalues[i] - eigenvalues[i - 1]) > 0.0001
            
            
        # end
    # end
    # plot!(plt, thetas, fs_norms, color = :black, title = latexstring("Eigenmodes of the \$\\gamma\$ Band"), label = "Fermi Surface")
    # display(plt)

end

main()