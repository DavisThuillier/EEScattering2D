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

    ## Collinear width parameters ##
    const max_width::Float64 = 0.4 # Maximum width of collinear region in radians
    const min_width::Float64 = 0.2 
    const width_slope::Float64 = 1.0 # Initial slope of exponential plateau function
    #####################
    
    collinear_width(T::Float64) = max_width * (1 - (1 - min_width/max_width) * exp( - (T - 0.0025) * width_slope))

    function gaussian_density(x::Float64, sigma::Float64, amplitude::Float64, limit_ratio::Float64)
        res = amplitude * exp( - (x / sigma)^2)
        if res <= limit_ratio * amplitude
            return limit_ratio * amplitude
        else 
            return res
        end 
    end 

    "Amplitude for density function normalized for N points between 0 point and endpoint."
    function get_amp(N::Float64, sigma::Float64, ratio::Float64, endpoint::Float64)
        return N / (sqrt(pi) * sigma * erf(log(1/ratio)) / 2 + (endpoint - sigma * log(1/ratio)) * ratio)
    end

    function get_energy_root(startpoint::SVector{2, Float64}, endpoint::SVector{2, Float64}, hamiltonian::Function; tolerance::Float64 = 0.0001, max_iterations::Int = 10000, level::Float64 = 0.0)
        delta_E::Float64 = 0.0
        n = endpoint - startpoint
        n = n / norm(n)

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
        fermi_surface = Vector{SVector{2,Float64}}(undef, length(angles))
        fill_fermi_surface!(fermi_surface, angles, hamiltonian; bz = bz)

        s = get_arclengths(fermi_surface)
        step::Float64 = last(s) / num_points

        uniform_fermi_surface = Vector{SVector{2, Float64}}(undef, num_points)
        uniform_fermi_surface[1] = fermi_surface[1]

        for i in 2:num_points
            t = step * (i - 1)
            j::Int = 1
            while t > s[j]
                j += 1
            end

            k = fermi_surface[j - 1] + ( (t - s[j - 1]) / (s[j] - s[j-1])) * (fermi_surface[j] - fermi_surface[j - 1])
            uniform_fermi_surface[i] = k
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
        arclengths = Vector{Float64}(undef, length(curve))
        arclengths[1] = 0.0

        for i in eachindex(curve)
            i == 1 || (arclengths[i] = arclengths[i - 1] + norm(curve[i] - curve[i - 1]) )
        end

        return arclengths
    end

    "Perimeter assuming curve is open using a straight line between the first and last points for closure"
    function get_perimeter(curve::Vector{SVector{2, Float64}})
        perimeter::Float64 = 0.0
        num_points::Int = lastindex(curve)
        for i in eachindex(curve)
            i < num_points ? (perimeter += norm(curve[i+1] - curve[i])) : (perimeter += norm(first(curve) - last(curve)))
        end
        return perimeter
    end

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

    function get_angles(fermi_surface::Vector{SVector{2, Float64}}, injection_index::Int, num_angles::Int, sigma::Float64)
        centered_fs = circshift(fermi_surface, 1 - injection_index) # Make the injection angle the first element of the array
        centered_fs_r = circshift(reverse(centered_fs), 1)
        t = get_arclengths(centered_fs) # Compute arclengths from the central angle
        t_r = get_arclengths(centered_fs_r)
        perimeter = get_perimeter(fermi_surface)

        sigma = sigma * perimeter / (2*pi) # Scale width of the peak by total arclength compared to unit circle
        central_angle = mod(atan(centered_fs[1][2], centered_fs[1][1]), 2 *pi)
        angles::Vector{Float64} = [central_angle]
        angles_r::Vector{Float64} = []
        limit_ratio = 0.2

        amp = get_amp(num_angles / 4, sigma, limit_ratio, t[end] / 4) 
        
        s::Float64     = 0.0
        i::Int         = 1
        theta::Float64 = central_angle
        k::Vector{Float64} = [0.0,0.0]

        extend_domain = (central_angle + pi >= 2*pi)
        while (theta < pi/2 + central_angle)
            s += 1 / gaussian_density(s, sigma, amp, limit_ratio)
            while s > t[i]
                i += 1
            end
            k = centered_fs[i - 1] + ( (s - t[i - 1]) / (t[i] - t[i-1])) * (centered_fs[i] - centered_fs[i - 1])
            theta = mod(atan(k[2], k[1]), 2*pi)
            extend_domain && ((theta <= pi) && (theta += 2*pi)) 
            (theta < pi/2 + central_angle) && push!(angles, theta)
        end
        push!(angles, pi/2 + central_angle)

        s = 0.0
        i = 1
        theta = central_angle
        extend_domain = (central_angle - pi <= 0)
        while theta > (central_angle - pi/2) && i < length(t_r)
            s += 1 / gaussian_density(s, sigma, amp, limit_ratio)
            while s > t_r[i]
                i += 1
            end
            k = centered_fs_r[i - 1] + ( (s - t_r[i - 1]) / (t_r[i] - t_r[i-1])) * (centered_fs_r[i] - centered_fs_r[i - 1])
            theta = mod(atan(k[2], k[1]), 2*pi)
            extend_domain && ((theta > pi) && (theta -= 2*pi))

            theta > (central_angle - pi/2) && push!(angles_r, theta)
        end

        return vcat(angles, pi .+ reverse(angles_r), pi .+ angles, reverse(angles_r))
    end

    function get_k_bound(hamiltonian::Function, e_bound::Float64, fs_k::SVector{2, Float64}, velocity::SVector{2,Float64}, max_iterations = 10000, tolerance = 0.0001; bz::Bool = true)
        e_bound  == 0.0 && return fs_k

        n = velocity / norm(velocity)
        step = (e_bound / norm(velocity)) / 2

        i::Int = 0
        i_limit::Int = abs(div(1, step)) # Number of step corresponding to half-width of Brillouin zone

        endpoint = fs_k
        startpoint = fs_k
        
        if e_bound > 0
            while i < i_limit && hamiltonian(endpoint) < e_bound
                endpoint += step * n
                bz && (abs(endpoint[1]) > 0.5 || abs(endpoint[2]) > 0.5) && break
                i += 1
            end
        else
            while i < i_limit && hamiltonian(startpoint) > e_bound
                startpoint += step * n
                i += 1
            end
        end

        j::Int = 0
        while j < max_iterations
            midpoint = (startpoint + endpoint) / 2
            delta_E = hamiltonian(midpoint) - e_bound
            norm(endpoint - startpoint) / 2 < tolerance && break
            
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
        # β::Float64 = sqrt(1 - 4*precision)
        # energies = Vector{Float64}(undef, perp_num)
        # energies[1] = 0.0
        # for i in 2:div(perp_num,2)+1
        #     energies[i] = T * log(1 / (fd(energies[i - 1], T) - β/perp_num) - 1)
        #     energies[perp_num - i + 2] = - energies[i] 
        # end
        # circshift!(energies, div(perp_num,2))
        #
        # for i in eachindex(fermi_surface)
        #     n = fermi_velocity[i] / norm(fermi_velocity[i])
        #     for j in 1:perp_num
        #         #momenta[j, i] = fermi_surface[i] + n * energies[j] 
        #         momenta[j,i] = get_k_bound(hamiltonian, energies[j], fermi_surface[i], fermi_velocity[i]; bz = bz)
        #     end
        # end


        ## Uniform Mesh ##
        for i in eachindex(fermi_surface)
            p_min = get_k_bound(hamiltonian, -e_max, fermi_surface[i], fermi_velocity[i]; bz = bz)
            p_max = get_k_bound(hamiltonian, e_max, fermi_surface[i], fermi_velocity[i]; bz = bz)

            for j in 1:perp_num
                @inbounds momenta[j, i] = p_min + (p_max - p_min) * (j - 1) / (perp_num - 1)
            end
        end

        return momenta
    end

    area(v::SVector{2, Float64}, u::SVector{2, Float64}) = abs(v[1] * u[2] - v[2] * u[1]) # Area in the plane spanned by 2-vectors

    function get_dV(ind_t::Int, ind_s::Int, momenta::Matrix{SVector{2, Float64}}, modulus)
        vertex = momenta[ind_t, ind_s]

        dV::Float64 = 0.0

        s_neighbors = [momenta[ind_t, mod(ind_s, modulus) + 1], momenta[ind_t, mod(ind_s - 2, modulus) + 1]]
        t_neighbors = [momenta[ind_t + 1, ind_s], momenta[ind_t - 1, ind_s]]
        diags = Matrix{SVector{2,Float64}}(undef, 2, 2)
        diags[1,1] = momenta[ind_t + 1, mod(ind_s, modulus) + 1]
        diags[1,2] = momenta[ind_t - 1, mod(ind_s, modulus) + 1]
        diags[2,1] = momenta[ind_t + 1, mod(ind_s - 2, modulus) + 1]
        diags[2,2] = momenta[ind_t - 1, mod(ind_s - 2, modulus) + 1]

        for i in eachindex(s_neighbors) 
            for j in eachindex(t_neighbors)
                mid_s = (s_neighbors[i] - vertex) / 2
                mid_t = (t_neighbors[j] - vertex) / 2
                diag = (diags[i,j] - vertex) / 2
                Δ1 = diag + (mid_t - mid_s)
                Δ2 = diag - (mid_t - mid_s)

                k1 = Δ1[2] / Δ1[1]
                k2 = Δ2[2] / Δ2[1]

                b1 = - k1 * mid_s[1] + mid_s[2]
                b2 = - k2 * mid_t[1] + mid_t[2]

                d::SVector{2, Float64} = [b1 + b2, k1*k2*(mid_t[1] - mid_s[1]) + k2 * mid_s[2] - k1 * mid_t[2]]
                
                dV += area(d, mid_t - mid_s) / 2
            end
        end 

        return dV
        
    end

    function discretize(fermi_surface::Vector{SVector{2, Float64}}, num_bins::Int, perp_num::Int, injection_index::Int, hamiltonian::Function, T::Float64, precision::Float64 = 0.001; bz::Bool = true)
        colin_width::Float64 = collinear_width(T)
        
        angles = mod2pi.( get_angles(fermi_surface, injection_index, num_bins, colin_width / (2 * sqrt(2 * log(2)))) ) # Get new angles on which to discretize

        arclengths = get_arclengths(fermi_surface)

        new_fs = Vector{SVector{2, Float64}}(undef, length(angles))
        new_fv = Vector{SVector{2, Float64}}(undef, length(angles))
        fill_fermi_surface!(new_fs, angles, hamiltonian; bz = bz)
        fill_fermi_velocity!(new_fv, new_fs, hamiltonian)

        dVs     = zeros(Float64, perp_num, length(new_fs))
        momenta = temperature_broaden(new_fs, new_fv, hamiltonian, perp_num, T, precision; bz = bz)

        points = map(x -> Point(x[1], x[2]), vec(momenta))
        BrillouinZone = Rectangle(Point(-0.5, -0.5), Point(0.5, 0.5))
        tess = voronoicells(points, BrillouinZone)
        areas = voronoiarea(tess)

        dVs = reshape(areas, size(momenta))
        # fill!(dVs[begin,:], 0.0)
        # fill!(dVs[end,:], 0.0)

        # for i in 2:(size(momenta)[1] - 1) # Ignore boundary contours
        #     for j in 1:size(momenta)[2]
        #         #dVs[i,j] = area((momenta[i + 1, j] - momenta[i - 1, j])/ 2, (momenta[i, mod(j, length(new_fs)) + 1] - momenta[i, mod(j - 2, length(new_fs)) + 1]) / 2)
        #         #dVs[i,j] = norm(momenta[i + 1, j] - momenta[i - 1, j])  * (norm(momenta[i, mod(j, length(new_fs)) + 1] - momenta[i,j]) + norm(momenta[i, mod(j - 2, length(new_fs)) + 1] - momenta[i,j])) / 4
        #         dVs[i,j] = get_dV(i,j, momenta, length(new_fs))
        #     end
        # end

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

        return momenta, dVs, variance, arclengths[injection_index] .+ get_arclengths(new_fs)
    end

    #####################
    ### Old Functions ###
    #####################

    # Stores the indices in the momentum array associated to a given angle
    struct Bin
        angle::Float64 
        indices::Vector{Tuple{Int, Int}}
    end

    function discretize_old(fermi_surface::Vector{SVector{2, Float64}}, num_bins::Int, perp_num::Int, injection_index::Int, hamiltonian::Function, T::Float64, precision::Float64 = 0.001)
        colin_width::Float64 = collinear_width(T)
        angles = mod2pi.( get_angles(fermi_surface, injection_index, num_bins, colin_width / (2 * sqrt(2 * log(2)))) ) # Get new angles on which to discretize

        new_fs = Vector{SVector{2, Float64}}(undef, length(angles))
        new_fv = Vector{SVector{2, Float64}}(undef, length(angles))
        fill_fermi_surface!(new_fs, angles, hamiltonian)
        fill_fermi_velocity!(new_fv, new_fs, hamiltonian)

        dV      = zeros(Float64, perp_num, length(new_fs))
        momenta = temperature_broaden(new_fs, new_fv, hamiltonian, perp_num, T, precision)
        
        angle_bins = Vector{Bin}(undef, length(angles))
        for i in eachindex(angles)
            angle_bins[i] = Bin(angles[i], [])
        end

        for i in 2:(size(momenta)[1] - 1) # Ignore boundary contours
            for j in 1:size(momenta)[2]
                dV[i,j] = area((momenta[i + 1, j] - momenta[i - 1, j])/ 2, (momenta[i, mod(j, length(new_fs)) + 1] - momenta[i, mod(j - 2, length(new_fs)) + 1]) / 2)

                theta = mod2pi(atan(momenta[i,j][2], momenta[i,j][1])) 
                min_index = argmin( abs.( - theta .+ angles))
                push!(angle_bins[min_index].indices, (i,j))
            end
        end

        for bin in angle_bins
            sort!(bin.indices, by = (x) -> norm(momenta[x[1], x[2]]))
        end

        variance = median(dV) / 4 # Approximate minimal width squared for delta function normalization

        return momenta, dV, angle_bins, variance
    end
end