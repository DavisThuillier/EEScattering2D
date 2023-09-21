# [src/mesh.jl]
module FermiSurfaceMesh
    export fill_fermi_surface!, fill_fermi_velocity!, discretize

    import StaticArrays: SVector
    import LinearAlgebra: norm
    import Statistics: median, mean
    import SpecialFunctions: erfinv, erf
    using ProgressBars
    using VoronoiCells
    import GeometryBasics: Point
    import StatsBase: countmap

    ## Collinear width parameters ##
    const max_width::Float64 = 0.4 # Maximum width of collinear region in radians
    const min_width::Float64 = 0.2 
    const width_slope::Float64 = 1.0 # Initial slope of exponential plateau function
    #####################
    
    # collinear_width(T::Float64) = max_width * (1 - (1 - min_width/max_width) * exp( - (T - 0.0025) * width_slope))
    const zero_temp_width = 1.84112 / 600
    const peak_slope      = 47.2989 / 600
    collinear_width(T::Float64) = zero_temp_width + peak_slope * T

    function get_energy_root(startpoint::SVector{2, Float64}, endpoint::SVector{2, Float64}, hamiltonian::Function; tolerance::Float64 = 0.0001, max_iterations::Int = 10000, level::Float64 = 0.0)
        delta_E::Float64 = 0.0

        j::Int = 0
        while j < max_iterations
            midpoint = (startpoint + endpoint) / 2
            delta_E = hamiltonian(midpoint) - level
            if norm(endpoint - startpoint) / 2 < tolerance
                return midpoint
            end
            
            if sign(delta_E) == sign(hamiltonian(startpoint))
                startpoint = midpoint
            else
                endpoint   = midpoint
            end

            j += 1
        end

        return (startpoint + endpoint) / 2
    end

    function generate_fermi_surface(hamiltonian::Function, num_points::Int; bz::Bool = true)
        angles = collect(range(0.0, 2*pi, num_points * 10))
        pop!(angles)
        fermi_surface = Vector{SVector{2,Float64}}(undef, length(angles))
        fill_fermi_surface!(fermi_surface, angles, hamiltonian; bz = bz)

        s = get_arclengths(fermi_surface)
        step::Float64 = last(s) / num_points

        uniform_fermi_surface = Vector{SVector{2, Float64}}(undef, num_points)
        uniform_fermi_surface[1] = fermi_surface[1]

        for i in 2:num_points
            t = step * (i - 1)
            uniform_fermi_surface[i] = get_momentum(fermi_surface, s, t)
        end

        return uniform_fermi_surface
    end

    function fill_fermi_surface!(fermi_surface::Vector{SVector{2, Float64}}, angles::Vector{Float64}, hamiltonian::Function; bz = true)
        startpoint::SVector{2, Float64} = [0.0, 0.0]
        center_energy = hamiltonian(startpoint)

        ## Compute roots of the Hamiltonian using the bisection method ##
        for i in eachindex(angles)
            n = SVector{2}([cos(angles[i]), sin(angles[i])])
            if bz
                if 0.0 <= angles[i] < pi / 4 || 3pi / 4 < angles[i] < 5pi / 4 || 7pi / 4 < angles[i] <= 2pi
                    endpoint   = 0.5 * sqrt(1 + sin(angles[i])^2) * n
                else
                    endpoint   = 0.5 * sqrt(1 + cos(angles[i])^2) * n
                end
            else
                endpoint = n
                j::Int = 1
                while sign(hamiltonian(endpoint)) == sign(center_energy) && j < 10
                    endpoint = n * j 
                    j += 1
                end
            end
            
            fermi_surface[i] = get_energy_root(startpoint, endpoint, hamiltonian, level = 0.0)
        end

        return nothing
    end

    function gradient(f::Function, k::SVector{2,Float64})
        dp::Float64 = sqrt(eps(Float64))
        df_x = f(k + SVector{2}([dp,0])) - f(k + SVector{2}([-dp, 0]))
        df_y = f(k + SVector{2}([0,dp])) - f(k + SVector{2}([0, -dp]))
        return SVector{2}([df_x, df_y] / (2 * dp))
    end

    function fill_fermi_velocity!(fermi_velocity::Vector{SVector{2, Float64}}, fermi_surface::Vector{SVector{2, Float64}}, hamiltonian::Function)
        for i in eachindex(fermi_surface)
            fermi_velocity[i] = gradient(hamiltonian, fermi_surface[i])
        end
        return nothing
    end

    function get_arclengths(curve::Vector{SVector{2, Float64}})
        arclengths = Vector{Float64}(undef, length(curve) + 1)
        arclengths[1] = 0.0

        for i in eachindex(curve)
            i == 1 && continue
            arclengths[i] = arclengths[i - 1] + norm(curve[i] - curve[i - 1])
        end
        arclengths[end] = arclengths[end - 1] + norm(curve[begin] - curve[end])

        return arclengths
    end

    function get_perimeter(curve::Vector{SVector{2, Float64}})
        return last(get_arclengths(curve))
    end

    "Assuming the curve is closed."
    function get_ds(curve::Vector{SVector{2, Float64}})
        ds = Vector{Float64}(undef, length(curve))
        if length(curve) > 1
            for i in eachindex(curve)
                ds[i] = 0.5 * ( norm(curve[mod(i, length(curve)) + 1] - curve[i]) + norm(curve[i] - curve[mod(i-2, length(curve)) + 1]))
            end
        else
            println("Curve must consist of at least two points.")
            exit()
        end
        return ds
    end

    function secant_method(f::Function, x0::Float64, x1::Float64, iterations::Int, precision::Float64)
        x2::Float64 = 0.0
        for i in 1:iterations
            x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
            x0, x1 = x1, x2
            abs(x0 - x1) < precision && break
        end 
        return x2
    end

    function minimum_binary_search(list::Vector{Float64}, target::Float64)
        left = 1
        right = length(list)
        while left <= right
            mid = div(left + right, 2)
            if list[mid] < target
                left = mid + 1
            elseif list[mid] > target
                right = mid - 1
            else
                return mid, mid
            end
        end
        return left - 1, left
    end

    function get_momentum(fermi_surface::Vector{SVector{2,Float64}}, arclengths::Vector{Float64}, s::Float64)
        i, j = mod.(minimum_binary_search(arclengths, s) .- 1, length(fermi_surface)) .+ 1
        if i == j
            return fermi_surface[i]
        else
            return fermi_surface[i] + (fermi_surface[j] - fermi_surface[i]) * (s - arclengths[i]) / (arclengths[j] - arclengths[i])
        end
    end

    # cdf(x::Float64, locus::Float64, amp::Float64, width::Float64, ratio::Float64) = amp * ((x - locus) / (width * sqrt(pi) * (ratio - 1)) + 0.5 * erf((x - locus) / width))

    cdf(x::Float64, x0::Float64, locus::Float64, width::Float64, ratio::Float64) = ((x - x0) / (width * sqrt(pi) * (ratio - 1)) + 0.5 * ( erf((x - locus) / width) - erf((x0 - locus) / width)) )

    function endpoint_gaussian_mesh()

    end


    function get_gaussian_fs(fermi_surface::Vector{SVector{2, Float64}}, loci::Vector{Float64}, num_angles::Int, width::Float64)
        arclengths = get_arclengths(fermi_surface)
        perimeter = arclengths[end]

        loci_indices = zeros(Int, 2 * length(loci))
        loci::Vector{Float64} = mod.(vcat(loci, loci .+ perimeter/2), perimeter)
        permutation = sortperm(loci)
        loci = loci[permutation]

        walls = Vector{Float64}(undef, length(loci)) # Midpoints between loci
        for i in eachindex(loci)
            i == length(loci) && continue
            walls[i] = mod((loci[i + 1] + loci[i]) / 2, perimeter)
            loci[i + 1] > loci[i] && (walls[end] += perimeter / 2)
        end
        walls[end] = (loci[1] + perimeter + loci[end]) / 2
        
        regions = Vector{Tuple{Float64, Float64}}(undef, 2*length(loci))
        for i in eachindex(regions)
            if isodd(i)
                regions[i] = (loci[div(i, 2) + 1], walls[div(i, 2) + 1])
            else
                regions[i] = (walls[div(i,2)], loci[mod(div(i, 2), length(loci)) + 1])
            end
        end
        regions[end] = (regions[end][1], regions[end][2] + perimeter)

        s_points = Vector{Float64}(undef, 0)

        ratio = 8.0
        # A = ((num_angles - length(loci)) / perimeter) / (1 + (ratio * length(loci)) / (perimeter) * erf(perimeter / (2 * length(loci) * width)))
        A = (0.5 * num_angles/length(loci)) / cdf(0.5 * perimeter / length(loci), 0.0, 0.0, width, ratio)
        N = 100 # Secant Method iteration limit

        for i in eachindex(regions)
            isodd(i) && (loci_indices[div(i + 1, 2)] = lastindex(s_points) + 1)
            isapprox(regions[i][2], regions[i][1]) && continue # Skip region if loci overlap

            if i != length(regions)
                locus = loci[div(i,2) + 1]
                push!(s_points, regions[i][1])
            else
                locus = regions[i][2]
                push!(s_points, locus)
            end
            s_points[end] - s_points[1] >= perimeter / 2 && break
            
            # push!(s_points, regions[i][1])
            s = locus
            j = 1
            if isodd(i)
                while j < num_angles
                    s = secant_method(x -> A * cdf(x, regions[i][1], locus, width, ratio) - j, regions[i][2], locus, N, 1e-8)
                    s - s_points[1] > perimeter / 2 && break
                    s >= regions[i][2] && break
                    push!(s_points, s)
                    j += 1
                end 
            else
                while j < num_angles
                    s = secant_method(x -> A * cdf(x, regions[i][2], locus, width, ratio) + j, locus, regions[i][1], N, 1e-8)
                    s - s_points[1] > perimeter / 2 && break
                    regions[i][1] >= s && break 
                    push!(s_points, s)
                    j += 1
                end
            end
        end
        sort!(s_points)
        
        mesh = Vector{SVector{2,Float64}}(undef, length(s_points) - 1)
        for i in eachindex(mesh)
            @inbounds mesh[i] = get_momentum(fermi_surface, arclengths, mod(s_points[i], perimeter))
        end

        return mesh, s_points[1] #loci_indices[invperm(permutation)][1:div(length(loci_indices), 2)]
    end

    function get_k_bound(hamiltonian::Function, e_bound::Float64, fs_k::SVector{2, Float64}, velocity::SVector{2,Float64}; max_iterations = 10000, tolerance = 0.0001, bz::Bool = true)
        e_bound  == 0.0 && return fs_k

        n = velocity / norm(velocity)
        step = (e_bound / norm(velocity))

        i::Int = 0
        i_limit::Int = abs(div(0.5, step)) # Number of step corresponding to half-width of Brillouin zone

        endpoint = fs_k
        startpoint = fs_k
        
        if e_bound > 0
            while hamiltonian(endpoint) < e_bound && i < i_limit
                endpoint += step * n
                bz && (abs(endpoint[1]) > 0.5 || abs(endpoint[2]) > 0.5) && break
                i += 1
            end
            endpoint += step * n
        else
            while hamiltonian(startpoint) > e_bound && i < i_limit
                startpoint += step * n
                i += 1
            end
            startpoint += step * n 
        end

        j::Int = 0
        while j < max_iterations
            midpoint = (startpoint + endpoint) / 2
            delta_E = hamiltonian(midpoint) - e_bound
            norm(startpoint - endpoint) < tolerance && break
            
            if sign(delta_E) == sign(hamiltonian(startpoint) - e_bound)
                startpoint = midpoint
            else
                endpoint   = midpoint
            end
            j += 1
        end

        k_bound = (startpoint + endpoint) / 2

        # Check if k_bound lies outside the Brillouin zone
        if bz
            if abs(k_bound[1]) > 0.5
                k_bound = fs_k + ( (0.5 * sign(k_bound[1]) - fs_k[1]) / n[1]) * n
            elseif abs(k_bound[2]) > 0.5
                k_bound = fs_k + ( (0.5 * sign(k_bound[2]) - fs_k[2]) / n[2]) * n
            end
        end

        return k_bound
    end

    fd(E::Float64, T::Float64) = 1 / (exp(E/T) + 1)

    function temperature_broaden(fermi_surface::Vector{SVector{2, Float64}}, fermi_velocity::Vector{SVector{2, Float64}}, hamiltonian::Function, perp_num::Int, T::Float64, precision::Float64; bz::Bool = true)
        e_max::Float64 = 2  * T * acosh(1 / (2 * sqrt(precision)))

        momenta = Matrix{SVector{2, Float64}}(undef, perp_num, length(fermi_surface))

        ## Fermi Profile ##
        β::Float64 = sqrt(1 - 4*precision)
        energies = Vector{Float64}(undef, perp_num)
        energies[1] = 0.0
        for i in 2:div(perp_num,2)+1
            @inbounds energies[i] = T * log(1 / (fd(energies[i - 1], T) - β/(perp_num-1)) - 1)
            energies[perp_num - i + 2] = - energies[i] 
        end
        
        circshift!(energies, div(perp_num,2))
        
        tolerance = abs(energies[2] / e_max) * T * 1e-2
        for i in eachindex(fermi_surface)            
            for j in 1:perp_num
                @inbounds momenta[j,i] = get_k_bound(hamiltonian, energies[j], fermi_surface[i], fermi_velocity[i]; tolerance = tolerance, bz = bz)
            end
        end

        return momenta
    end

    function discretize(fermi_surface::Vector{SVector{2, Float64}}, num_bins::Int, perp_num::Int, loci::Vector{Float64}, hamiltonian::Function, T::Float64, precision::Float64 = 0.001; bz::Bool = true)
        perimeter = get_perimeter(fermi_surface)
        width::Float64 = 4.0 * collinear_width(T) / (2 * sqrt(log(2))) * perimeter # Really, σ sqrt(2) for the Gaussian density with a 5 σ minimal peak width

        half_fs, s_primus = get_gaussian_fs(fermi_surface, loci, num_bins, width)
        half_fv = Vector{SVector{2, Float64}}(undef, length(half_fs))
        fill_fermi_velocity!(half_fv, half_fs, hamiltonian)
        
        momenta = temperature_broaden(half_fs, half_fv, hamiltonian, perp_num, T, precision; bz = bz)
        momenta = hcat(momenta, - momenta) # Enforce inversion symmetry

        # @show map(x -> findall(isequal(x), momenta), collect(keys(filter(kv -> kv.second > 1, countmap(momenta)))) )

        arclengths = get_arclengths(momenta[div(perp_num,2) + 1, :])
        arclengths = s_primus .+ arclengths

        size(unique(momenta)) != size(vec(momenta)) && (println("Duplicate points generated in mesh."); return nothing)

        dVs     = zeros(Float64, size(momenta))

        points = map(x -> Point(x[1], x[2]), vec(momenta))
        BrillouinZone = Rectangle(Point(-0.5, -0.5), Point(0.5, 0.5))
        tess = voronoicells(points, BrillouinZone)
        areas = voronoiarea(tess)

        dVs = reshape(areas, size(momenta))
        fill!(dVs[begin,:], 0.0)
        fill!(dVs[end,:], 0.0)

        mean_dE = 0.0
        for j in 1:size(momenta)[2]
            min_dE = 1.0 # Set to be the width of the BZ 
            dE = 0.0
            for i in 1:(size(momenta)[1] - 1)
                dE = abs(hamiltonian(momenta[i + 1,j]) - hamiltonian(momenta[i, j]))
                dE < min_dE && (min_dE = dE)
            end
            mean_dE += min_dE
        end
        mean_dE = mean_dE / size(momenta)[2]

        variance = mean_dE^2 / 2 # Approximate minimal width squared for delta function normalization

        loci_indices = Vector{Int}(undef, length(loci))
        for i in eachindex(loci)
            j1, j2 = minimum_binary_search(arclengths, mod(loci[i], perimeter))
            loci_indices[i] = div(j1 + j2, 2)
        end

        return momenta, dVs, variance, arclengths, loci_indices
    end

    discretize(fermi_surface::Vector{SVector{2, Float64}}, num_bins::Int, perp_num::Int, locus::Float64, hamiltonian::Function, T::Float64, precision::Float64 = 0.001; bz = true) = discretize(fermi_surface, num_bins, perp_num, [locus], hamiltonian, T, precision; bz = bz)
    
    discretize(fermi_surface::Vector{SVector{2, Float64}}, num_bins::Int, perp_num::Int, injection_index::Int, hamiltonian::Function, T::Float64, precision::Float64 = 0.001; bz = true) = discretize(fermi_surface, num_bins, perp_num, [get_arclengths(fermi_surface)[injection_index]], hamiltonian, T, precision; bz = bz)

    function endpoint_gaussian_mesh(a::Real, b::Real, N::Int, amplitude_ratio::Real, width::Real)
        b < a && ((a, b) = (b, a))
    
        loci = [a, (a + b) / 2.0, b]
        
        amp = (N/ 4.0 - 1) / cdf((b - a) / 4.0, 0.0, 0.0, width, amplitude_ratio)
        
        mesh = Vector{Float64}(undef, 0)
        domain = LinRange(0.0, (b - a) / 2.0, 10 * Int(N * amplitude_ratio))
    
        for i in 1:2
            if isodd(i)
                integrated_n = amp * cdf.(domain, first(domain), first(domain), width, amplitude_ratio)
            else
                integrated_n = amp * cdf.(domain, first(domain), last(domain), width, amplitude_ratio)
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
    
        return vcat(mesh, [b])
    end

end