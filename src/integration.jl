# [src/integration.jl]
module FermiSurfaceIntegration
    export collision_integral, contracted_integral!

    import StaticArrays: SVector
    import LinearAlgebra: norm

    import ..EEScattering2D: FermiSurfaceMesh

    fd(E::Float64, T::Float64) = 1 / (exp(E/T) + 1)

    gaussian_delta(deviation::Float64, width::Float64) = exp( - deviation^2 / width) / sqrt(width * pi)

    function gradient(f::Function, k::SVector{2,Float64})
        dp::Float64 = sqrt(eps(Float64))
        df_x = f(k + SVector{2}([dp,0])) - f(k + SVector{2}([-dp, 0]))
        df_y = f(k + SVector{2}([0,dp])) - f(k + SVector{2}([0, -dp]))
        return SVector{2}([df_x, df_y] / (2 * dp))
    end

    function collision_integral(p1::SVector{2,Float64}, k::SVector{2,Float64}, integration_mesh::Matrix{SVector{2, Float64}}, energies::Matrix{Float64}, dVs::Matrix{Float64}, hamiltonian::Function, sigma_squared::Float64, T::Float64, q_squared::Float64; umklapp = true)
        
        I2_n::Float64 = 0.0
        I2_u::Float64 = 0.0
        I34_n::Float64 = 0.0
        I34_u::Float64 = 0.0

        P = p1 + k # Total momentum for process 2
        Q = p1 - k # Momentum difference for p1 -> p1'
        
        E1 = hamiltonian(p1)
        Ek = hamiltonian(k)
        E_sum::Float64 = E1 + Ek
        E_diff::Float64 = E1 - Ek

        mod_shift = SVector{2}([0.5, 0.5]) # For shifting the wavevector over before taking the modulus with respect to the first Brillouin Zone

        for i in 2:(size(integration_mesh)[1] - 1)
            for j in 1:size(integration_mesh)[2]
                @inbounds p1_prime = integration_mesh[i,j]
                p2_prime = umklapp ? mod.(mod_shift + P - p1_prime, 1.0) - mod_shift : P - p1_prime
                E1_prime = energies[i,j]
                E2_prime = hamiltonian(p2_prime)

                # V = 0.5 * (1 / (norm(p1_prime - p1)^2 + q_squared) + 1 / (norm(p1_prime - k)^2 + q_squared))
                if abs((P - p1_prime)[1]) > 0.5 || abs((P - p1_prime)[2]) > 0.5 # Detect whether this is an umklapp event
                    I2_u += fd(Ek, T) * (1 - fd(E1_prime, T)) * (1 - fd(E2_prime, T)) * dVs[i, j] * gaussian_delta(E_sum - E1_prime - E2_prime, sigma_squared) #/ (norm(p1_prime - p1)^2 + q_squared)^2
                else
                    I2_n += fd(Ek, T) * (1 - fd(E1_prime, T)) * (1 - fd(E2_prime, T)) * dVs[i, j] * gaussian_delta(E_sum - E1_prime - E2_prime, sigma_squared) # / (norm(p1_prime - p1)^2 + q_squared)^2
                end

                @inbounds p2 = integration_mesh[i,j]
                p2_prime = umklapp ? mod.(mod_shift + Q + p2, 1.0) - mod_shift : Q + p2
                E2 = energies[i,j]
                E2_prime = hamiltonian(p2_prime)
                if abs((Q + p2)[1]) > 0.5 || abs((Q + p2)[2]) > 0.5 # Detect whether this is an umklapp event
                    I34_u += fd(E2, T) * (1 - fd(Ek, T)) * (1 - fd(E2_prime, T)) * dVs[i,j] * gaussian_delta((E_diff + E2 - E2_prime), sigma_squared) #/ (norm(p1 - k)^2 + q_squared)^2
                else
                    I34_n += fd(E2, T) * (1 - fd(Ek, T)) * (1 - fd(E2_prime, T)) * dVs[i,j] * gaussian_delta((E_diff + E2 - E2_prime), sigma_squared)
                end

                @inbounds p2_prime = integration_mesh[i,j]
                p2 = umklapp ? mod.(mod_shift + p2_prime - Q, 1.0) - mod_shift : p2_prime - Q
                E2_prime = energies[i,j]
                E2 = hamiltonian(p2)
                if abs((Q + p2)[1]) > 0.5 || abs((Q + p2)[2]) > 0.5 # Detect whether this is an umklapp event
                    I34_u += fd(E2, T) * (1 - fd(Ek, T)) * (1 - fd(E2_prime, T)) * dVs[i,j] * gaussian_delta((E_diff + E2 - E2_prime), sigma_squared) #/ (norm(p1 - k)^2 + q_squared)^2
                else
                    I34_n += fd(E2, T) * (1 - fd(Ek, T)) * (1 - fd(E2_prime, T)) * dVs[i,j] * gaussian_delta((E_diff + E2 - E2_prime), sigma_squared) #/ (norm(p1 - k)^2 + q_squared)^2
                end

            end
        end
        
        return [(- I2_n + I34_n) * fd(E1, T),  (- I2_u + I34_u) * fd(E1, T)] # fd(E1, T) is a constant in each integral and thus removed
    end

    # function collision_integral_old(p1_index::Tuple{Int, Int}, k_index::Tuple{Int,Int}, integration_mesh::Matrix{SVector{2, Float64}}, energies::Matrix{Float64}, dVs::Matrix{Float64}, hamiltonian::Function, sigma_squared::Float64, T::Float64, q_squared::Float64; umklapp = true)
    #     p1::SVector{2,Float64} = integration_mesh[p1_index[1], p1_index[2]]
    #     k::SVector{2,Float64} = integration_mesh[k_index[1], k_index[2]]
        
    #     I2_n::Float64 = 0.0
    #     I2_u::Float64 = 0.0
    #     I34_n::Float64 = 0.0
    #     I34_u::Float64 = 0.0

    #     P = p1 + k # Total momentum of scatterers
    #     Q = p1 - k # Momentum difference for p1 -> p1'

    #     E1 = hamiltonian(p1)
    #     Ek = hamiltonian(k)
    #     E_sum::Float64 = E1 + Ek
    #     E_diff::Float64 = E1 - Ek

    #     mod_shift = SVector{2}([0.5, 0.5]) # For shifting the wavevector over before taking the modulus with respect to the first Brillouin Zone

    #     for i in 2:(size(integration_mesh)[1] - 1)
    #         for j in 1:size(integration_mesh)[2]
    #             @inbounds p1_prime = integration_mesh[i,j]
    #             p2_prime = umklapp ? mod.(mod_shift + P - p1_prime, 1.0) - mod_shift : P - p1_prime
    #             E1_prime = energies[i,j]
    #             E2_prime = hamiltonian(p2_prime)
    #             if abs((P - p1_prime)[1]) > 0.5 || abs((P - p1_prime)[2]) > 0.5 # Detect whether this is an umklapp event
    #                 I2_u += fd(Ek, T) * (1 - fd(E1_prime, T)) * (1 - fd(E2_prime, T)) * dVs[i, j] * gaussian_delta(E_sum - E1_prime - E2_prime, sigma_squared) / (norm(p1_prime - p1)^2 + q_squared)^2
    #             else
    #                 I2_n += fd(Ek, T) * (1 - fd(E1_prime, T)) * (1 - fd(E2_prime, T)) * dVs[i, j] * gaussian_delta(E_sum - E1_prime - E2_prime, sigma_squared) / (norm(p1_prime - p1)^2 + q_squared)^2
    #             end

    #             @inbounds p2 = integration_mesh[i,j]
    #             p2_prime = umklapp ? mod.(mod_shift + Q + p2, 1.0) - mod_shift : Q + p2
    #             E2 = energies[i,j]
    #             E2_prime = hamiltonian(p2_prime)
    #             if abs((Q + p2)[1]) > 0.5 || abs((Q + p2)[2]) > 0.5 # Detect whether this is an umklapp event
    #                 I34_u += 2 * fd(E2, T) * (1 - fd(Ek, T)) * (1 - fd(E2_prime, T)) * dVs[i,j] * gaussian_delta((E_diff + E2 - E2_prime), sigma_squared) / (norm(p1 - k)^2 + q_squared)^2
    #             else
    #                 I34_n += 2 * fd(E2, T) * (1 - fd(Ek, T)) * (1 - fd(E2_prime, T)) * dVs[i,j] * gaussian_delta((E_diff + E2 - E2_prime), sigma_squared) / (norm(p1 - k)^2 + q_squared)^2
    #             end
    #         end
    #     end
        
    #     return [(- I2_n + I34_n) * fd(E1, T),  (- I2_u + I34_u) * fd(E1, T)] # fd(E1, T) is a constant in each integral and thus removed
    # end

    "Compute Boltzmann collision integral between each point on FS and the injection point."

    function contracted_integral!(reduced_mat::Matrix{Float64}, arclengths::Vector{Float64}, perimeter::Float64, uniform_fs::Vector{SVector{2,Float64}}, momenta::Matrix{SVector{2, Float64}}, hamiltonian::Function, T::Float64, q_squared::Float64; umklapp = true)
        integral::SVector{2,Float64} = [0.0,0.0]
        s_num = size(momenta)[2]
        t_num = size(momenta)[1]

        central_momenta::Vector{SVector{2, Float64}} = momenta[:, 1]
        central_dp = zeros(t_num)

        for i in 2:(t_num - 1)
            central_dp[i] = norm(central_momenta[i + 1] - central_momenta[i - 1]) / 2
        end

        for i in 1:s_num            
            integral = SVector{2}([0.0, 0.0])
            integration_mesh, mesh_dVs, mesh_variance, _ , _ = FermiSurfaceMesh.discretize(uniform_fs, s_num, t_num, [arclengths[begin], arclengths[i]], hamiltonian, T, 0.001)
            energies = hamiltonian.(integration_mesh) 

            for m in 2:(t_num-1)
                k = momenta[m, i]
                dp = norm(momenta[m + 1, i] - momenta[m - 1, i]) / 2
                for j in eachindex(central_momenta)
                    p1 = momenta[j, 1]
                    loss_terms = collision_integral(p1, k, integration_mesh, energies, mesh_dVs, hamiltonian, mesh_variance, T, q_squared; umklapp = umklapp) / ((2pi)^4 * T)

                    integral += loss_terms * dp * central_dp[j] 
                end                  
            end
            reduced_mat[i + 1, 1] = mod.(arclengths[i], perimeter)
            reduced_mat[i + 1, 2] = integral[1] #/ width # Normal scattering contribution
            reduced_mat[i + 1, 3] = integral[2] #/ width # Umklapp scattering contribution
        end

        ### Padding rows ###
        reduced_mat[begin, :] = [ -3*perimeter, 0.0, 0.0] # Will be sorted as the first element
        reduced_mat[end, :] = [ 3*perimeter, 0.0, 0.0] # Will be sorted as the final element
        return nothing
    end
end