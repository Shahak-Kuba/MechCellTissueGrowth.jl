"""
    Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt, type)

Calculate the normal velocity of the interface such that `Vₙ` is proportional to `ρ`. The direction of the normal vector is determined by `type`.

# Arguments
- `rᵢ₋₁`: The point before the current point in space.
- `rᵢ`: The current point in space.
- `rᵢ₊₁`: The point after the current point in space.
- `kf`: The amount of tissue produced per unit area per unit time.
- `δt`: The time step.
- `type`: A string specifying the orientation of the normal vector, either "inward" or any other value for outward orientation.

# Returns
The normal velocity of the interface.
"""
function Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf, δt,type)
    ρₗ = ρ(rᵢ, rᵢ₋₁)
    ρᵣ = ρ(rᵢ₊₁, rᵢ)
    Vₗ = kf .* ρₗ
    Vᵣ = kf .* ρᵣ

    nₗ = n(rᵢ₋₁, rᵢ,type)
    nᵣ = n(rᵢ, rᵢ₊₁,type)

    rₘ₁ = rᵢ - (Vₗ .* nₗ .* δt)
    rₗ = rᵢ₋₁ - (Vₗ .* nₗ .* δt)
    rₘ₂ = rᵢ - (Vᵣ .* nᵣ .* δt)
    rᵣ = rᵢ₊₁ - (Vᵣ .* nᵣ .* δt)

    return (lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ) - rᵢ) ./ δt
end


function Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, ρ, n, kf, δt)
    ρₗ = circshift(ρ,(0,1))
    Vₗ = (kf .* ρₗ')'
    Vᵣ = (kf .* ρ')'

    nₗ = circshift(n,(0,1))

    rₘ₁ = rᵢ .+ (Vₗ .* nₗ .* δt)
    rₗ = rᵢ₋₁ .+ (Vₗ .* nₗ .* δt)
    rₘ₂ = rᵢ .+ (Vᵣ .* n .* δt)
    rᵣ = rᵢ₊₁ .+ (Vᵣ .* n .* δt)

    return (lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ) .- rᵢ) ./ δt
end



# alternate normal velocity function with no corrections based on spring overlaps or seperations
function Vₙ(rᵢ₋₁, rᵢ, rᵢ₊₁, kf,type)
    ρₗ = ρ(rᵢ, rᵢ₋₁)
    ρᵣ = ρ(rᵢ₊₁, rᵢ)
    V = kf .* ((ρₗ .+ ρᵣ)./2)

    nᵥ = n(rᵢ₊₁,rᵢ₋₁,type)

    return V.*nᵥ
end


