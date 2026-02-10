# User defined force function

## Hookean Restoring force
hookean_restoring_force = (l, kₛ, a) -> (kₛ .* (l' .- a))'
## Nonlinear restoring force
nonlinear_restoring_force = (l, kₛ, a) -> (kₛ .* (1 ./ a .- 1 ./ l'))'
nonlinear_restoring_force2 = (l, kₛ, a) -> (kₛ .* a.^2 .* (1 ./ a .- 1 ./ l'))'


# When changing force law make sure to run all of these
#FORCE_FNC = (rᵢ, rⱼ, kₛ, l₀) -> nonlinear_restoring_force(rᵢ, rⱼ, kₛ, l₀)

# Force functions used in ODEs
#Fₛ⁺(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) =  FORCE_FNC(rᵢ, rᵢ₊₁, kₛ, l₀) .* τ(rᵢ₊₁, rᵢ)
#Fₛ⁻(rᵢ, rᵢ₊₁, rᵢ₋₁, kₛ, l₀) = -FORCE_FNC(rᵢ, rᵢ₋₁, kₛ, l₀) .* τ(rᵢ, rᵢ₋₁)

function F(l, kₛ, a, restoring_force)
    if restoring_force == "hookean"
        return (kₛ .* (l' .- a))'
    elseif restoring_force == "nonlinear"
        return (kₛ .* (1 ./ a .- 1 ./ l'))'
    else
        return (kₛ .* a.^2 .* (1 ./ a .- 1 ./ l'))'
    end
end


function Fₛ⁺(l⁺, kₛ, a, restoring_force)
    force_function = get_force_function(restoring_force)
    return force_function(l⁺, kₛ, a)
end

function Fₛ⁻(l⁻, kₛ, a, restoring_force)
    force_function = get_force_function(restoring_force)
    return -force_function(l⁻, kₛ, a)
end

function get_force_function(restoring_force)
    if restoring_force == "hookean"
        return hookean_restoring_force
    elseif restoring_force == "nonlinear"
        return nonlinear_restoring_force
    else
        return nonlinear_restoring_force2
    end
end

