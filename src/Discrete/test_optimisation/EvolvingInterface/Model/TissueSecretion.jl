"""
    Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt, type)

Calculate normal velocity using line-intersection method (original interface).
Used by post-processing code and the old-style ODE; the optimised ODE inlines this.
"""
function Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt, type)
    ρₗ = ρ(rᵢ, rᵢ₋₁)
    ρᵣ = ρ(rᵢ₊₁, rᵢ)
    Vₗ = kf .* ρₗ
    Vᵣ = kf .* ρᵣ

    nₗ = n(rᵢ₋₁, rᵢ, type)
    nᵣ = n(rᵢ, rᵢ₊₁, type)

    rₘ₁ = rᵢ .- (Vₗ .* nₗ .* δt)
    rₗ  = rᵢ₋₁ .- (Vₗ .* nₗ .* δt)
    rₘ₂ = rᵢ .- (Vᵣ .* nᵣ .* δt)
    rᵣ  = rᵢ₊₁ .- (Vᵣ .* nᵣ .* δt)

    return (lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ) .- rᵢ) ./ δt
end

"""
    Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, ρ, n, kf, δt)

Vectorised normal velocity using pre-computed densities and normals.
Used by PostCalcs2D.
"""
function Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, ρ, n, kf, δt)
    ρₗ = circshift(ρ, (0, 1))
    Vₗ = (kf .* ρₗ')'
    Vᵣ = (kf .* ρ')'

    nₗ = circshift(n, (0, 1))

    rₘ₁ = rᵢ .+ (Vₗ .* nₗ .* δt)
    rₗ  = rᵢ₋₁ .+ (Vₗ .* nₗ .* δt)
    rₘ₂ = rᵢ .+ (Vᵣ .* n .* δt)
    rᵣ  = rᵢ₊₁ .+ (Vᵣ .* n .* δt)

    return (lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ) .- rᵢ) ./ δt
end

"""Simple normal velocity without line-intersection correction."""
function Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, type)
    ρₗ = ρ(rᵢ, rᵢ₋₁)
    ρᵣ = ρ(rᵢ₊₁, rᵢ)
    V  = kf .* ((ρₗ .+ ρᵣ) ./ 2)
    nᵥ = n(rᵢ₊₁, rᵢ₋₁, type)
    return V .* nᵥ
end
