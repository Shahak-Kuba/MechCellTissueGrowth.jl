# =====================================================================
#  Force functions
# =====================================================================

# Vectorised force lambdas (used by post-processing / analysis)
hookean_restoring_force     = (l, kₛ, a) -> (kₛ .* (l' .- a))'
nonlinear_restoring_force   = (l, kₛ, a) -> (kₛ .* (1 ./ a .- 1 ./ l'))'
nonlinear_restoring_force2  = (l, kₛ, a) -> (kₛ .* a.^2 .* (1 ./ a .- 1 ./ l'))'

"""Vectorised force computation dispatched on `ForceType` enum (no string comparison)."""
function F(l, kₛ, a, ft::ForceType)
    if ft == Hookean
        return (kₛ .* (l' .- a))'
    elseif ft == Nonlinear
        return (kₛ .* (1 ./ a .- 1 ./ l'))'
    else  # Nonlinear2
        return (kₛ .* a.^2 .* (1 ./ a .- 1 ./ l'))'
    end
end

# Backward-compat wrapper for FB_ODE! which still passes a string
function F(l, kₛ, a, ft::String)
    ftype = ft == "hookean" ? Hookean : (ft == "nonlinear" ? Nonlinear : Nonlinear2)
    return F(l, kₛ, a, ftype)
end

"""Positive spring force (from current node towards next node)."""
function Fₛ⁺(l⁺, kₛ, a, ft::ForceType)
    if ft == Hookean
        return hookean_restoring_force(l⁺, kₛ, a)
    elseif ft == Nonlinear
        return nonlinear_restoring_force(l⁺, kₛ, a)
    else
        return nonlinear_restoring_force2(l⁺, kₛ, a)
    end
end

"""Negative spring force (from current node towards previous node)."""
function Fₛ⁻(l⁻, kₛ, a, ft::ForceType)
    if ft == Hookean
        return -hookean_restoring_force(l⁻, kₛ, a)
    elseif ft == Nonlinear
        return -nonlinear_restoring_force(l⁻, kₛ, a)
    else
        return -nonlinear_restoring_force2(l⁻, kₛ, a)
    end
end
