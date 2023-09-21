

module NewMesh

    import StaticArrays: SVector 
    import SpecialFunctions: erf
    import Statistics: median, mean
    using VoronoiCells
    import GeometryBasics: Point

    ## Collinear width parameters ##
    const max_width::Float64 = 0.4 # Maximum width of collinear region in radians
    const min_width::Float64 = 0.2 
    const width_slope::Float64 = 1.0 # Initial slope of exponential plateau function
    #####################
    
    collinear_width(T::Float64) = max_width * (1 - (1 - min_width/max_width) * exp( - (T - 0.0025) * width_slope))

    function get_arclengths(curve::Vector{SVector{2, Float64}})
        arclengths = Vector{Float64}(undef, length(curve) + 1)
        arclengths[1] = 0.0

        for i in eachindex(curve)
            i == 1 && continue
            arclengths[i] = arclengths[i - 1] + norm(curve[i] - curve[i - 1])
        end
        arclengths[end] = arclengths[end - 1] + norm(curve[begin] - curve[end]) # Close the curve to get the perimeter

        return arclengths
    end

    function get_perimeter(curve::Vector{SVector{2, Float64}})
        return last(get_arclengths(curve))
    end

    function integrated_gaussian_density(x::Real, x0::Real, amp_ratio::Real, width::Real)
        return 0.5 * (erf((x - x0) / width) - erf(-x0 / width)) + x / (width * sqrt(pi) * (amp_ratio - 1.0))
    end

    function endpoint_gaussian_mesh(a::Real, b::Real, N::Int, amplitude_ratio::Real, width::Real)
        b < a && ((a, b) = (b, a))
    
        loci = [a, (a + b) / 2.0, b]
        
        amp = (N/ 4.0 - 1) / integrated_gaussian_density((b - a) / 4.0, 0.0, amplitude_ratio, width)
        
        mesh = Vector{Float64}(undef, 0)
        domain = LinRange(0.0, (b - a) / 2.0, 100 * round(Int, N * amplitude_ratio))
    
        for i in 1:2
            if isodd(i)
                integrated_n = amp * integrated_gaussian_density.(domain, first(domain),amplitude_ratio, width)
            else
                integrated_n = amp * integrated_gaussian_density.(domain, last(domain),amplitude_ratio, width)
            end
            n::Int = 0
            push!(mesh, loci[i])
            for j in eachindex(integrated_n)
                j == 1 && continue
                if integrated_n[j] > n
                    push!(mesh, loci[i] + domain[j] + (n - integrated_n[j - 1]) * (domain[j] - domain[j - 1]) / (integrated_n[j] - integrated_n[j - 1])) 
                    n += 1
                end
            end
        end
        push!(mesh, b)
    
        return mesh
        
    end

    function fs_gaussian_mesh(perimeter::Real, loci::Vector{<:Real}, N::Int, amplitude_ratio::Real, width::Real)
        loci_indices = zeros(Int, 2 * length(loci))
        loci::Vector{Float64} = mod.(vcat(loci, loci .+ perimeter/2), perimeter)
        permutation = sortperm(loci)
        loci = loci[permutation]

        @show loci

        mesh = Vector{Real}(undef, 0)
        for i in eachindex(loci)
            i == 1 && continue
            region_points = round(Int, N * (loci[i] - loci[i - 1]) / perimeter)
            mesh = vcat(mesh, endpoint_gaussian_mesh(loci[i-1], loci[i], region_points, amplitude_ratio, width))
        end

        return mesh
    end
    
    function gradient(f::Function, k::SVector{2,Float64})
        dp::Float64 = sqrt(eps(Float64))
        df_x = f(k + SVector{2}([dp,0])) - f(k + SVector{2}([-dp, 0]))
        df_y = f(k + SVector{2}([0,dp])) - f(k + SVector{2}([0, -dp]))
        return SVector{2}([df_x, df_y] / (2 * dp))
    end

    fd(E::Float64, T::Float64) = 1 / (exp(E/T) + 1)

    # function temperature_broaden(fermi_surface::Vector{SVector{2, Float64}}, fermi_velocity::Vector{SVector{2, Float64}}, hamiltonian::Function, perp_num::Int, T::Float64, precision::Float64)
    #     e_max::Float64 = 2  * T * acosh(1 / (2 * sqrt(precision)))

    #     momenta = Matrix{SVector{2, Float64}}(undef, perp_num, length(fermi_surface))

    #     ## Fermi Profile ##
    #     β::Float64 = sqrt(1 - 4*precision)
    #     energies = Vector{Float64}(undef, perp_num)
    #     energies[1] = 0.0
    #     for i in 2:div(perp_num,2)+1
    #         @inbounds energies[i] = T * log(1 / (fd(energies[i - 1], T) - β/(perp_num-1)) - 1)
    #         energies[perp_num - i + 2] = - energies[i] 
    #     end
        
    #     circshift!(energies, div(perp_num,2))
        
    #     tolerance = abs(energies[2] / e_max) * T * 1e-2
    #     for i in eachindex(fermi_surface)            
    #         for j in 1:perp_num
    #             @inbounds momenta[j,i] = get_k_bound(hamiltonian, energies[j], fermi_surface[i], fermi_velocity[i]; tolerance = tolerance)
    #         end
    #     end

    #     return momenta
    # end

    # function discretize(fermi_surface::Vector{SVector{2, Float64}}, num_bins::Int, perp_num::Int, loci::Vector{Float64}, hamiltonian::Function, T::Float64, precision::Float64 = 0.001)
    #     width::Float64 = collinear_width(T) / (2 * sqrt(log(2))) # Really, σ sqrt(2) for the Gaussian density

    #     arclengths = get_arclengths(fermi_surface)

        
    #     loci::Vector{Float64} = mod.(vcat(loci, loci .+ perimeter/2), perimeter)


    #     # perimeter = get_perimeter(fermi_surface)

        
    #     # half_fs, s_primus = get_gaussian_fs(fermi_surface, loci, num_bins, sigma)
    #     # half_fv = Vector{SVector{2, Float64}}(undef, length(half_fs))
    #     # fill_fermi_velocity!(half_fv, half_fs, hamiltonian)
        
    #     # momenta = temperature_broaden(half_fs, half_fv, hamiltonian, perp_num, T, precision)
    #     # momenta = hcat(momenta, - momenta) # Enforce inversion symmetry

    #     # # @show map(x -> findall(isequal(x), momenta), collect(keys(filter(kv -> kv.second > 1, countmap(momenta)))) )

    #     # arclengths = get_arclengths(momenta[div(perp_num,2) + 1, :])
    #     # arclengths = s_primus .+ arclengths

    #     # size(unique(momenta)) != size(vec(momenta)) && (println("Duplicate points generated in mesh."); return nothing)

    #     # dVs     = zeros(Float64, size(momenta))

    #     # points = map(x -> Point(x[1], x[2]), vec(momenta))
    #     # BrillouinZone = Rectangle(Point(-0.5, -0.5), Point(0.5, 0.5)) # In units of (2pi / a)
    #     # tess = voronoicells(points, BrillouinZone)
    #     # areas = voronoiarea(tess)

    #     # dVs = reshape(areas, size(momenta))
    #     # fill!(dVs[begin,:], 0.0)
    #     # fill!(dVs[end,:], 0.0)

    #     # mean_dE = 0.0
    #     # for j in 1:size(momenta)[2]
    #     #     min_dE = 1.0 # Set to be the width of the BZ 
    #     #     dE = 0.0
    #     #     for i in 1:(size(momenta)[1] - 1)
    #     #         dE = abs(hamiltonian(momenta[i + 1,j]) - hamiltonian(momenta[i, j]))
    #     #         dE < min_dE && (min_dE = dE)
    #     #     end
    #     #     mean_dE += min_dE
    #     # end
    #     # mean_dE = mean_dE / size(momenta)[2]

    #     # variance = mean_dE^2 / 2 # Approximate minimal width squared for delta function normalization

    #     # loci_indices = Vector{Int}(undef, length(loci))
    #     # for i in eachindex(loci)
    #     #     j1, j2 = minimum_binary_search(arclengths, mod(loci[i], perimeter))
    #     #     loci_indices[i] = div(j1 + j2, 2)
    #     # end

    #     return nothing
    # end

end