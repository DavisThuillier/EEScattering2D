# [src/integration.jl]
module FermiSurfaceIntegration
    export collision_integral, angular_integral!

    import StaticArrays: SVector
    import LinearAlgebra: norm

    fd(E::Float64, T::Float64) = 1 / (exp(E/T) + 1)

    fd_normalization(E::Float64, T::Float64) = 4 * cosh(E / (2 * T))

    gaussian_delta(deviation::Float64, sigma_squared::Float64) = exp( - deviation^2 / (2 * sigma_squared)) / (sigma_squared * sqrt(2 * pi))

    function collision_integral(p1::SVector{2, Float64}, k::SVector{2, Float64}, momenta::Matrix{SVector{2, Float64}}, dVs::Matrix{Float64}, hamiltonian::Function, sigma_squared::Float64, T::Float64)
        I2_n::Float64 = 0.0 
        I2_u::Float64 = 0.0
        I34_n::Float64 = 0.0
        I34_u::Float64 = 0.0

        P = p1 + k # Total momentum of scatterers
        Q = p1 - k # Momentum difference for p1 -> p1'

        E2::Float64 = hamiltonian(p1) + hamiltonian(k)
        E34::Float64 = hamiltonian(p1) - hamiltonian(k)

        mod_shift = SVector{2}([1.0,1.0]) # For shifting the wavevector over before taking the modulus with respect to the first Brillouin Zone

        for i in 2:(size(momenta)[1] - 1)
            for j in 1:size(momenta)[2]
                @inbounds p1_prime = momenta[i,j]
                p2_prime = mod.(mod_shift + P - p1_prime, 2.0) - mod_shift
                #p2_prime = P - p1_prime
                if ( abs((P - p1_prime)[1] > 1) || abs((P - p1_prime)[1] > 1) ) # Detect whether this is an umklapp event
                    I2_u += fd(hamiltonian(k), T) * (1 - fd(hamiltonian(p1_prime), T)) * (1 - fd(hamiltonian(p2_prime), T)) * dVs[i] * gaussian_delta(E2 - hamiltonian(p1_prime) - hamiltonian(p2_prime), sigma_squared) 
                else
                    I2_n += fd(hamiltonian(k), T) * (1 - fd(hamiltonian(p1_prime), T)) * (1 - fd(hamiltonian(p2_prime), T)) * dVs[i] * gaussian_delta(E2 - hamiltonian(p1_prime) - hamiltonian(p2_prime), sigma_squared) 
                end

                @inbounds p2 = momenta[i,j]
                p2_prime = mod.(mod_shift + Q + p2, 2) - mod_shift
                #p2_prime = Q + p2
                if ( abs((Q + p2)[1] > 1) || abs((Q + p2)[1] > 1) ) # Detect whether this is an umklapp event
                    I34_n += 2 * fd(hamiltonian(p2), T) * (1 - fd(hamiltonian(k), T)) * (1 - fd(hamiltonian(p2_prime), T)) * dVs[i] * gaussian_delta((E34 + hamiltonian(p2) - hamiltonian(p2_prime)), sigma_squared)
                else
                    I34_u += 2 * fd(hamiltonian(p2), T) * (1 - fd(hamiltonian(k), T)) * (1 - fd(hamiltonian(p2_prime), T)) * dVs[i] * gaussian_delta((E34 + hamiltonian(p2) - hamiltonian(p2_prime)), sigma_squared)
                end
            end
        end
        
        return [(- I2_n + I34_n) * fd(hamiltonian(p1), T),  (- I2_u + I34_u) * fd(hamiltonian(p1), T)] # fd(p1, T) is a constant in each integral and thus removed
    end

    "Compute Boltzmann collision integral between each point on FS and the injection point."
    function contracted_integral!(reduced_mat::Matrix{Float64}, arclengths::Vector{Float64}, perimeter::Float64, momenta::Matrix{SVector{2, Float64}}, dVs::Matrix{Float64}, hamiltonian::Function, variance::Float64, T::Float64)
        integral::SVector{2,Float64} = [0.0,0.0]

        perp_num = size(momenta)[1]
        central_momenta::Vector{SVector{2, Float64}} = momenta[:, 1]
        central_dp = zeros(perp_num)

        for i in 2:(perp_num - 1)
            central_dp[i] = norm(central_momenta[i + 1] - central_momenta[i - 1]) / 2
        end
          
        for i in 1:size(momenta)[2]            
            integral = SVector{2}([0.0, 0.0])
            for m in 2:(perp_num-1)
                k = momenta[m, i]
                dp = norm(momenta[m + 1, i] - momenta[m - 1, i]) / 2
                for j in eachindex(central_momenta)
                    p1 = central_momenta[j]
                    loss_terms = collision_integral(p1, k, momenta, dVs, hamiltonian, variance, T) * fd_normalization(hamiltonian(p1), T)

                    integral += loss_terms * dp * central_dp[j] / T # Divide by T due to delta functions in the integration
                end   
            end
            reduced_mat[i + 1, 1] = mod.(arclengths[i], perimeter)
            reduced_mat[i + 1, 2] = integral[1] # Normal scattering contribution
            reduced_mat[i + 1, 3] = integral[2] # Umklapp scattering contribution
        end

        ### Padding rows ###
        reduced_mat[begin, :] = [ - 3*perimeter, 0.0, 0.0] # Will be sorted as the first element
        reduced_mat[end, :] = [ 3*perimeter, 0.0, 0.0] # Will be sorted as the final element
        return nothing
    end

    #####################
    ### Old Functions ###
    #####################

    # function get_dp(momenta::Matrix{SVector{2,Float64}}, bin::Bin)
    #     # Returns the radial momentum differential assuming sorted momenta in the bin
    #     bin_momenta = map((x) -> momenta[x[1], x[2]], bin.indices)
    #     dp = Vector{Float64}(undef, sizeof(bin_momenta))

    #     for i in eachindex(bin_momenta)
    #         if i < length(bin_momenta)
    #             dp[i] = norm(bin_momenta[i + 1]) - norm(bin_momenta[i])
    #         else 
    #             dp[i] = dp[i - 1]           
    #         end
    #     end

    #     return dp
    # end

    # "Compute Boltzmann collision integral between each binned angle and the injection angle."
    # function angular_integral!(reduced_mat::Matrix{Float64}, momenta::Matrix{SVector{2, Float64}}, dVs::Matrix{Float64}, bins::Vector{Bin}, hamiltonian::Function, variance::Float64, T::Float64)
    #     integral::SVector{2,Float64} = [0.0,0.0]

    #     central_momenta::Vector{SVector{2, Float64}} = map((x) -> momenta[x[1], x[2]], bins[1].indices)
    #     central_dp::Vector{Float64} = get_dp(momenta, bins[1])
        
    #     for i in eachindex(bins)
    #         dp = get_dp(momenta, bins[i])
    #         integral = SVector{2}([0.0, 0.0])
    #         for m in eachindex(bins[i].indices)
    #             coord = bins[i].indices[m]
    #             k = momenta[coord[1], coord[2]]
    #             for j in eachindex(central_momenta)
    #                 p1 = central_momenta[j]
    #                 loss_terms = collision_integral(p1, k, momenta, dVs, hamiltonian, variance, T) * fd_normalization(hamiltonian(p1), T)

    #                 integral += loss_terms * (norm(p1) * norm(k) * dp[m] * central_dp[j]) / T # Divide by T due to delta functions in the integration
    #             end   
    #         end
    #         reduced_mat[i, 1] = bins[i].angle
    #         reduced_mat[i, 2] = integral[1] # Normal scattering contribution
    #         reduced_mat[i, 3] = integral[2] # Umklapp scattering contribution
    #     end

    #     ### Padding rows ###
    #     reduced_mat[end - 1, :] = [-6*pi, 0.0, 0.0] # Will be sorted as the first element
    #     reduced_mat[end, :] = [6*pi, 0.0, 0.0] # Will be sorted as the final element
    #     return nothing
    # end

end