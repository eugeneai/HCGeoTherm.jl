# Empirical Geotherm Calculations Module
# Based on method by Hasterok & Chapman [EPSL, 2011]
# Refactored to remove global variables and improve code structure

"""
    empgtherms(q0, maxsz, dz, D, zbot, H)

Compute empirical geotherms based on Hasterok & Chapman [EPSL, 2011] method.

# Arguments
- `q0`: Surface heat flow [mW/m²]
- `maxsz`: Maximum depth of model [km]
- `dz`: Depth step (layer width) [km]
- `D`: Thickness of upper crust [km]
- `zbot`: Depth to layer base [km] (vector)
- `H`: Layer heat production [μW/m³] (vector)

# Returns
- `T`: Temperature as a function of depth [°C]
- `z`: Depth [km]
- `k`: Thermal conductivity [W/m/K]
- `A`: Heat production [W/m³]
- `q`: Heat flow [W/m²]
- `alpha`: Thermal expansivity [1/K]
- `az`: Adiabat depth if reached, else -1.0 [km]

# Example
```julia
T, z, k, A, q, alpha, az = empgtherms(40.0, 225.0, 0.1, 16.0, [16, 23, 39, 300], [0.0, 0.4, 0.4, 0.02])
```
"""
function empgtherms(q0::Float64, maxsz::Float64, dz::Float64, D::Float64,
                    zbot::Vector{Float64}, H::Vector{Float64})

    # Validate inputs
    @assert length(zbot) == length(H) "zbot and H must have same length"
    @assert all(zbot .> 0) "zbot must be positive"
    @assert all(H .>= 0) "H must be non-negative"
    @assert q0 > 0 "q0 must be positive"
    @assert maxsz > 0 "maxsz must be positive"
    @assert dz > 0 "dz must be positive"
    @assert D > 0 "D must be positive"

    # Depths
    z = 0:dz:maxsz
    lenz = length(z)

    # Physical constants
    g = 10.0  # gravity [m/s²]
    rhoc = 2850.0  # crust density [kg/m³]
    rhom = 3340.0  # mantle density [kg/m³]
    zmoho = zbot[3]  # Moho depth [km]

    # Initialize arrays
    lambda = zeros(lenz - 1)
    alpha = zeros(lenz)

    # Compute heat production profile
    A = compute_heat_production_profile(z, dz, D, zbot, H)

    # Temperature conditions
    T0 = 0.0 + 273.15  # Surface temperature [K]
    Ta0 = 1300.0 + 273.15  # Adiabatic temperature at surface [K]
    dT = 0.3  # Adiabatic gradient [K/km]

    # Adiabatic gradient
    Ta = Ta0 .+ z .* dT

    # Compute heat flow
    q = compute_heat_flow(q0, A, dz, lenz)

    # Compute temperature profile
    T, lambda, alpha, adiabat_reached, adiabat_depth, last_index =
        compute_temperature_profile(z, dz, T0, Ta, q, lambda, alpha, A, zmoho, zbot, rhoc, rhom, g)

    # Handle adiabat if reached
    if adiabat_reached
        T, q, A, alpha, lambda = compute_adiabat_profile(
            T, Ta, q, A, alpha, lambda, z, last_index, lenz, zmoho, zbot, rhoc, rhom, g
        )
    end

    # Convert temperature to Celsius
    T = T .- 273.15

    return T, collect(z), lambda, A, q, alpha, adiabat_depth
end

"""
    compute_heat_production_profile(z, dz, D, zbot, H)

Compute heat production profile based on layer boundaries.
"""
function compute_heat_production_profile(z::Vector{Float64}, dz::Float64, D::Float64,
                                         zbot::Vector{Float64}, H::Vector{Float64})

    lenz = length(z)
    A = zeros(lenz)

    current_layer = 1
    layer_end = round(Int64, D / dz)

    for i in 1:length(zbot)
        if i != 1
            layer_start = layer_end + 1
            layer_end = round(Int64, zbot[i] / dz)
        else
            layer_start = 1
            layer_end = round(Int64, D / dz)
        end

        # Ensure indices are within bounds
        layer_start = min(layer_start, lenz)
        layer_end = min(layer_end, lenz)

        if layer_start <= layer_end
            A[layer_start:layer_end] .= H[i]
        end
    end

    return A
end

"""
    compute_heat_flow(q0, A, dz, lenz)

Compute heat flow profile from surface heat flow and heat production.
"""
function compute_heat_flow(q0::Float64, A::Vector{Float64}, dz::Float64, lenz::Int)
    q = fill(q0, lenz)

    # Heat flow decreases with depth due to heat production
    if lenz > 1
        cumulative_A = cumsum(A[1:lenz-1]) * dz
        q[2:lenz] = q[2:lenz] .- cumulative_A
    end

    return q
end

"""
    compute_temperature_profile(z, dz, T0, Ta, q, lambda, alpha, A, zmoho, zbot, rhoc, rhom, g)

Compute temperature profile until adiabat is reached or maximum depth.
"""
function compute_temperature_profile(z::Vector{Float64}, dz::Float64, T0::Float64, Ta::Vector{Float64},
                                    q::Vector{Float64}, lambda::Vector{Float64}, alpha::Vector{Float64},
                                    A::Vector{Float64}, zmoho::Float64, zbot::Vector{Float64},
                                    rhoc::Float64, rhom::Float64, g::Float64)

    lenz = length(z)
    T = zeros(lenz)
    T[1] = T0

    current_layer = 1
    adiabat_reached = false
    adiabat_depth = -1.0
    last_index = lenz

    # Initial thermal expansivity at surface
    alpha[1] = empexpansivity(current_layer, zmoho, z[1], T[1])

    for i in 1:lenz-1
        # Update current layer if depth exceeds layer boundary
        if z[i] >= zbot[current_layer] && current_layer < length(zbot)
            current_layer += 1
        end

        # Compute temperature and thermal conductivity
        T[i+1], lambda[i] = tccomp(
            current_layer, zmoho, z[i+1], dz, T[i], q[i],
            i == 1 ? 3.0 : lambda[i-1], A[i]
        )

        # Compute thermal expansivity
        alpha[i+1] = empexpansivity(current_layer, zmoho, z[i+1], T[i+1])

        # Check if adiabat is reached
        if T[i+1] > Ta[i+1]
            adiabat_reached = true
            adiabat_depth = z[i]
            last_index = i
            break
        end
    end

    return T, lambda, alpha, adiabat_reached, adiabat_depth, last_index
end

"""
    compute_adiabat_profile(T, Ta, q, A, alpha, lambda, z, start_idx, lenz, zmoho, zbot, rhoc, rhom, g)

Compute adiabatic profile from where geotherm reaches adiabat to maximum depth.
"""
function compute_adiabat_profile(T::Vector{Float64}, Ta::Vector{Float64}, q::Vector{Float64},
                                 A::Vector{Float64}, alpha::Vector{Float64}, lambda::Vector{Float64},
                                 z::Vector{Float64}, start_idx::Int, lenz::Int, zmoho::Float64,
                                 zbot::Vector{Float64}, rhoc::Float64, rhom::Float64, g::Float64)

    current_layer = find_layer(z[start_idx], zbot)

    for j in (start_idx+1):lenz
        T[j] = Ta[j]
        q[j] = q[j-1]
        A[j] = A[j-1]

        # Update current layer if depth exceeds layer boundary
        if z[j] >= zbot[current_layer] && current_layer < length(zbot)
            current_layer += 1
        end

        # Compute thermal expansivity
        alpha[j] = empexpansivity(current_layer, zmoho, z[j], T[j])

        if j == lenz
            break
        end

        # Compute thermal conductivity for adiabatic region
        lambda[j], _ = thermcond(current_layer, zmoho, z[j], 0.5 * (Ta[j] + Ta[j+1]))
    end

    return T, q, A, alpha, lambda
end

"""
    find_layer(z_depth, zbot)

Find the layer index for a given depth.
"""
function find_layer(z_depth::Float64, zbot::Vector{Float64})
    for (i, boundary) in enumerate(zbot)
        if z_depth < boundary
            return i
        end
    end
    return length(zbot)
end

"""
    tccomp(ik, zmoho, z, dz, tau, q, lambda0, A)

Compute temperature and thermal conductivity for a depth step using Newton-Raphson iteration.
"""
function tccomp(ik::Int, zmoho::Float64, z::Float64, dz::Float64,
                tau::Float64, q::Float64, lambda0::Float64, A::Float64)

    # Starting temperature guess
    dT0 = q / lambda0 * dz
    guess = tau + dT0

    max_iterations = 20
    tolerance = 1e-6
    gamma = q * dz - 0.5 * A * dz^2

    for iteration in 1:max_iterations
        # Thermal conductivity and its derivative
        lambda, dlambda = thermcond(ik, zmoho, z, 0.5 * (guess + tau))

        # Newton-Raphson iteration
        f = guess - (tau + gamma / lambda)
        df = 1 + gamma * dlambda / lambda^2

        new_guess = guess - f / df

        # Check convergence
        if abs(new_guess - guess) < tolerance
            guess = new_guess
            break
        end

        guess = new_guess

        if iteration == max_iterations
            @warn "TCCOMP did not converge after $max_iterations iterations"
        end
    end

    # Final thermal conductivity
    lambda_final, _ = thermcond(ik, zmoho, z, 0.5 * (guess + tau))

    return guess, lambda_final
end

"""
    thermcond(ik, zmoho, z, T)

Compute thermal conductivity and its temperature derivative.

# Arguments
- `ik`: Layer index
- `zmoho`: Moho depth [km]
- `z`: Depth [km]
- `T`: Temperature [K]

# Returns
- `lambda`: Thermal conductivity [W/m/K]
- `dlambda`: Derivative of thermal conductivity with respect to temperature
"""
function thermcond(ik::Int, zmoho::Float64, z::Float64, T::Float64)
    # Get thermal conductivity coefficients
    k = kcoef(ik, T)

    # Physical constants
    g = 10.0
    rhoc = 2850.0
    rhom = 3340.0

    # Compute pressure
    if z < zmoho
        P = 1e-6 * rhoc * g * z
    else
        P = 1e-6 * (rhoc * g * zmoho + rhom * g * (z - zmoho))
    end

    # Compute thermal conductivity and its derivative
    lambda = (k[1] + k[2] / T + k[3] * T^2) * (1 + k[4] * P)
    dlambda = (2 * k[3] * T - k[2] / T^2) * (1 + k[4] * P)

    return lambda, dlambda
end

"""
    kcoef(ik, T)

Get thermal conductivity coefficients for a given layer and temperature.

# Returns
Array of coefficients [a, b, c, d] where:
lambda = (a + b/T + c*T²) * (1 + d*P)
"""
function kcoef(ik::Int, T::Float64)
    # Coefficients for different rock types
    # Format: [a, b, c, d]

    # Upper crust coefficients (T < 844 K)
    ka = [
        [1.496,  398.84,  4.573e-7, 0.0950],  # Layer 1
        [1.733,  194.59,  2.906e-7, 0.0788],  # Layer 2
        [1.723,  219.88,  1.705e-7, 0.0520],  # Layer 3
        [2.271,  681.12, -1.259e-7, 0.0399],  # Layer 4
        [2.371,  669.40, -1.288e-7, 0.0384]   # Layer 5
    ]

    # Lower crust coefficients (T > 844 K)
    kb = [
        [2.964, -495.29,  0.866e-7, 0.0692],  # Layer 1
        [2.717, -398.93,  0.032e-7, 0.0652],  # Layer 2
        [2.320,  -96.98, -0.981e-7, 0.0463]   # Layer 3
    ]

    # Select coefficients based on temperature and layer
    if ik <= 3 && T > 844.0
        return kb[ik]
    else
        # Ensure ik is within bounds
        idx = min(ik, length(ka))
        return ka[idx]
    end
end

"""
    empexpansivity(ia, zmoho, z, T)

Compute thermal expansivity.

# Arguments
- `ia`: Layer index
- `zmoho`: Moho depth [km]
- `z`: Depth [km]
- `T`: Temperature [K]

# Returns
- `alpha`: Thermal expansivity [1/K]
"""
function empexpansivity(ia::Int, zmoho::Float64, z::Float64, T::Float64)
    # Get thermal expansivity coefficients
    a = acoef(ia, T)

    # Physical constants
    g = 10.0
    rhoc = 2850.0
    rhom = 3340.0

    # Compute pressure
    if z < zmoho
        P = 1e-6 * rhoc * g * z
    else
        P = 1e-6 * (rhoc * g * zmoho + rhom * g * (z - zmoho))
    end

    # Compute thermal expansivity
    alpha = (a[1] + a[2] * T + a[3] / T^2) * (1 + a[4] * P)

    return alpha
end

"""
    acoef(ia, T)

Get thermal expansivity coefficients for a given layer and temperature.

# Returns
Array of coefficients [a, b, c, d] where:
alpha = (a + b*T + c/T²) * (1 + d*P)
"""
function acoef(ia::Int, T::Float64)
    # Coefficients for different rock types
    # Format: [a, b, c, d]

    # Upper crust coefficients (T < 844 K)
    aa = [
        [2.355e-5, 3.208e-8, -0.7938, -0.1193],  # Layer 1
        [2.020e-5, 2.149e-8, -0.6315, -0.1059],  # Layer 2
        [2.198e-5, 0.921e-8, -0.1820, -0.0626],  # Layer 3
        [3.036e-5, 0.925e-8, -0.2730, -0.0421],  # Layer 4
        [3.026e-5, 0.906e-8, -0.3116, -0.0408]   # Layer 5
    ]

    # Lower crust coefficients (T > 844 K)
    ab = [
        [1.741e-5, 0.500e-8, -0.3094, -0.0778],  # Layer 1
        [1.663e-5, 0.602e-8, -0.3364, -0.0745],  # Layer 2
        [2.134e-5, 0.711e-8, -0.1177, -0.0563]   # Layer 3
    ]

    # Select coefficients based on temperature and layer
    if ia <= 3 && T > 844.0
        return ab[ia]
    else
        # Ensure ia is within bounds
        idx = min(ia, length(aa))
        return aa[idx]
    end
end

# Export functions
export empgtherms, tccomp, thermcond, kcoef, empexpansivity, acoef

end # module
