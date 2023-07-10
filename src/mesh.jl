# [src/mesh.jl]
module FermiSurfaceMesh

    export fill_fermi_surface!, fill_fermi_velocity!, discretize

    import StaticArrays: SVector

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

    function get_energy_root(startpoint::SVector{2, Float64}, endpoint::SVector{2, Float64}, hamiltonian::Function, level::Float64 = 0.0)
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

    function fill_fermi_surface!(fermi_surface::Vector{SVector{2, Float64}}, angles::Vector{Float64})
        startpoint::SVector{2, Float64} = [0.0, 0.0]

        ## Compute roots of the Hamiltonian using the bisection method ##
        for i in eachindex(angles)
            n = SVector{2}([cos(angles[i]), sin(angles[i])])
            endpoint = 2.0 * n 
            if - pi / 4 < angles[i] < pi / 4 || 3 * pi / 4 < angles[i] < 5 * pi / 4
                endpoint   = sqrt(1 + sin(angles[i])^2) * n
            else
                endpoint   = sqrt(1 + cos(angles[i])^2) * n
            end
            fermi_surface[i] = get_energy_root(startpoint, endpoint, hamiltonian, 0.0)
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
        arc_lengths = Vector{Float64}(undef, length(curve))
        arc_lengths[1] = 0.0

        for i in 2:length(curve)
            arc_lengths[i] = arc_lengths[i - 1] + norm(curve[i] - curve[i - 1])
        end

        return arc_lengths
    end

    function get_angles(fermi_surface::Vector{SVector{2, Float64}}, injection_index::Int, num_angles::Int, sigma::Float64)
        centered_fs = circshift(fermi_surface, 1 - injection_index) # Make the injection angle the first element of the array
        centered_fs_r = circshift(reverse(centered_fs), 1)
        t = get_arclengths(centered_fs) # Compute arc_lengths from the central angle
        t_r = get_arclengths(centered_fs_r)
        sigma = sigma * t[end] / (2*pi) # Scale width of the peak by total arclength compared to unit circle
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
        while (theta < pi/2 + central_angle) && i < length(t)
            s += 1 / gaussian_density(s, sigma, amp, limit_ratio)
            while s > t[i]
                i += 1
            end
            k = centered_fs[i - 1] + ( (s - t[i - 1]) / (t[i] - t[i-1])) * (centered_fs[i] - centered_fs[i - 1])
            theta = mod(atan(k[2], k[1]), 2*pi)
            extend_domain && ((theta <= pi) && (theta += 2*pi)) 
            (theta < pi/2 + central_angle) && push!(angles, theta)
        end

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

        (length(angles_r) > num_angles / 2) && err("Discretization error: too many angles.")
        (length(angles) > num_angles / 2) && err("Discretization error: too many angles.")

        return vcat(angles, pi .+ reverse(angles_r), pi .+ angles, reverse(angles_r))
    end

    function get_k_bound(e_bound::Float64, fs_k::SVector{2, Float64}, velocity::SVector{2,Float64})
        n = velocity / norm(velocity)
        step = (e_bound / norm(velocity)) / 2

        i::Int = 0
        i_limit::Int = abs(div(1, step)) # Number of step corresponding to half-width of Brillouin zone

        endpoint = fs_k
        startpoint = fs_k
        
        if e_bound > 0
            while i < i_limit && hamiltonian(endpoint) < e_bound
                endpoint += step * n
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

        #Check if k_bound lies outside the Brillouin zone
        if abs(k_bound[1]) > 1
            k_bound = fs_k + ( (sign(k_bound[1]) - fs_k[1]) / n[1]) * n
        elseif abs(k_bound[2]) > 1
            k_bound = fs_k + ( (sign(k_bound[2]) - fs_k[2]) / n[2]) * n
        end
        return k_bound
    end

    function temperature_broaden(fermi_surface::Vector{SVector{2, Float64}}, fermi_velocity::Vector{SVector{2, Float64}}, T::Float64)
        e_max::Float64 = 2  * T * acosh(1 / (2 * sqrt(prec)))

        momenta = Matrix{SVector{2, Float64}}(undef, p_num, length(fermi_surface))

        for i in eachindex(fermi_surface)
            # p_max = get_k_bound(e_max, n, fermi_surface[i], step)
            p_min = get_k_bound(-e_max, fermi_surface[i], fermi_velocity[i])
            p_max = get_k_bound(e_max, fermi_surface[i], fermi_velocity[i])

            for j in 1:p_num # 2 additional points correspond to boundary contours of ring
                @inbounds momenta[j, i] = p_min + (p_max - p_min) * (j - 1) / (p_num - 1)
            end
        end

        return momenta
    end

    area(v::SVector{2, Float64}, u::SVector{2, Float64}) = abs(v[1] * u[2] - v[2] * u[1]) # Area in the plane spanned by 2-vectors

    # Stores the indices in the momentum array associated to a given angle
    struct Bin
        angle::Float64 
        indices::Vector{Tuple{Int, Int}}
    end

    function discretize(fermi_surface::Vector{SVector{2, Float64}}, injection_index::Int, T::Float64)
        colin_width::Float64 = collinear_width(T)
        angles = mod2pi.( get_angles(fermi_surface, injection_index, num_angles, colin_width / (2 * sqrt(2 * log(2)))) ) # Get new angles on which to discretize

        new_fs = Vector{SVector{2, Float64}}(undef, length(angles))
        new_fv = Vector{SVector{2, Float64}}(undef, length(angles))
        fill_fermi_surface!(new_fs, angles)
        fill_fermi_velocity!(new_fv, new_fs, hamiltonian)

        dV      = zeros(Float64, p_num, length(new_fs))
        momenta = temperature_broaden(new_fs, new_fv, T)
        
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