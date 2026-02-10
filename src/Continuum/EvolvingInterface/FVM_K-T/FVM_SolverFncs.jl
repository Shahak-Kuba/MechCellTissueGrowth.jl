#using Polynomials


function minmod(varᵢ₋₁, varᵢ, varᵢ₊₁, Φ, Δx)
    # flux limiter function given by [Kraganov-Tadmor 2000] FVM Scheme equation (4.9), 
    # this function is a vectorised code based of the implementation seen in Alias & Buenzli 2017.
    s1 = sign.(Φ.*(varᵢ₊₁ .- varᵢ)./Δx)
    s2 = sign.((varᵢ₊₁ .- varᵢ₋₁)/(2 .*Δx))
    s3 = sign.(Φ.*(varᵢ .- varᵢ₋₁)./Δx)

    output = zeros(size(varᵢ))

    for i in eachindex(output)
        if s1[i] == s2[i] && s1[i] == s3[i]
            output[i] = s1[i]*min(abs(Φ*(varᵢ₊₁[i] - varᵢ[i])/Δx), abs((varᵢ₊₁[i] - varᵢ₋₁[i])/(2 *Δx)), abs(Φ*(varᵢ[i] - varᵢ₋₁[i])/Δx))
        end
    end

    return output
end


# Functions from Modified conservation laws
function f₁(r, σ, η)
    return zeros(size(r))
end

function f₂(r, σ, η, kf)
    return kf.*η./r
end

function f₃(r, σ, η, kf, D, dη, κ)
    return (kf.*σ.*(η.^2))./(r.^3 + r.*σ.^2) .- (D.*dη)./(r.^2 .+ σ.^2) .+ (D.*η.*κ.*σ)./(r.*sqrt.(r.^2 .+ σ.^2)) .+ (2 .*D.*η.*σ)./(r.^3 + r.*σ.^2)
end


# Functions to calculate eigenvalues of the Jacobian
function β(r, σ, η, kf, D, dη, dσ, ddσ)
    return (kf.*(η.^2).*r)./((r.^2 .+ σ.^2).^2) .- (kf.*(η.^2).*(σ.^2))./(r.*(r.^2 .+ σ.^2).^2) .+
            (D.*η.*r)./((r.^2 .+ σ.^2).^2) .- (8 .*D.*η.*(σ.^2))./(r.*(r.^2 .+ σ.^2).^2) .+
            (2 .*D.*σ.*dη)./((r.^2 .+ σ.^2).^2) .+ (4 .*D.*η.*(σ.^2).*r)./((r.^2 .+ σ.^2).^3) .+
            (8 .*D.*η.*(σ.^4))./(r.*(r.^2 .+ σ.^2).^3) .+ (D.*η.*dσ)./((r.^2 .+ σ.^2).^2) .-
            (4 .*D.*η.*(σ.^2).*dσ)./((r.^2 .+ σ.^2).^3) #.+ (D.*η.*σ.*ddσ)./((r.^2 .+ σ.^2).^2)
end

function Γ(r, σ, η, kf, D, dσ)
    return (2 .*kf.*σ.*η .+ 2 .*D.*σ)./(r.*(r.^2 .+ σ.^2)) .- 
           (D.*σ.*r)./((r.^2 .+ σ.^2).^2) .- 
           (2 .*D.*(σ.^3))./(r.*(r.^2 .+ σ.^2).^2) .+
           (D.*σ.*dσ)./((r.^2 .+ σ.^2).^2)
end

# Solving and returning max eigenvalues (i.e. roots to [-1 γ β/r 0] and characteristic equation given by "-λ³ + γλ² + (β/r)λ = 0" )
function λ(Β,Γ,r,kf)
    λ_max = zeros(size(r));

    λ₊ = (abs.((Γ .+ sqrt.(Complex.(Γ.^2 .- (4 .*kf.*Β)./r)))./-2))
    λ₋ = (abs.((Γ .- sqrt.(Complex.(Γ.^2 .- (4 .*kf.*Β)./r)))./-2))

    for i in eachindex(λ_max)
        λ_max[i] = max(λ₊[i],λ₋[i],0.0)
    end
    """ Method using Roots function (creates a lot of allocations)
    p_coefficients = hcat(zeros(length(λ_max)), kf .* Β ./ r, Γ, -ones(length(λ_max)))

    p_polynomials = Polynomial.(eachcol(p_coefficients'))
    λs = roots.(p_polynomials)
    

    for i in eachindex(λ_max)
        λ_max[i] = Float64(maximum(abs.(λs[i]))[1])
    end"""
    return λ_max
end

function a⁺(λ⁺₁, λ⁺₂)
    return max.(λ⁺₁, λ⁺₂)
end

function a⁻(λ⁻₁, λ⁻₂)
    return max.(λ⁻₁, λ⁻₂)
end

function SetupSolverMemory(r)
    # vectors used in simulation (for speed)
    r⁺₁ = zeros(size(r[1,:]))
    r⁺₂ = zeros(size(r[1,:]))
    σ⁺₁ = zeros(size(r[1,:]))
    σ⁺₂ = zeros(size(r[1,:]))
    η⁺₁ = zeros(size(r[1,:]))
    η⁺₂ = zeros(size(r[1,:]))
    r⁻₁ = zeros(size(r[1,:]))
    r⁻₂ = zeros(size(r[1,:]))
    σ⁻₁ = zeros(size(r[1,:]))
    σ⁻₂ = zeros(size(r[1,:]))
    η⁻₁ = zeros(size(r[1,:]))
    η⁻₂ = zeros(size(r[1,:]))
    dσ⁺₁ = zeros(size(r[1,:]))
    dσ⁺₂ = zeros(size(r[1,:]))
    dη⁺₁ = zeros(size(r[1,:]))
    dη⁺₂ = zeros(size(r[1,:]))
    dσ⁻₁ = zeros(size(r[1,:]))
    dσ⁻₂ = zeros(size(r[1,:]))
    dη⁻₁ = zeros(size(r[1,:]))
    dη⁻₂ = zeros(size(r[1,:]))
    β⁺₁ = zeros(size(r[1,:]))
    γ⁺₁ = zeros(size(r[1,:]))
    β⁺₂ = zeros(size(r[1,:]))
    γ⁺₂ = zeros(size(r[1,:]))
    λ⁺₁ = zeros(size(r[1,:]))
    λ⁺₂ = zeros(size(r[1,:]))
    a⁺val = zeros(size(r[1,:]))
    β⁻₁ = zeros(size(r[1,:]))
    γ⁻₁ = zeros(size(r[1,:]))
    β⁻₂ = zeros(size(r[1,:]))
    γ⁻₂ =zeros(size(r[1,:]))
    λ⁻₁ = zeros(size(r[1,:]))
    λ⁻₂ = zeros(size(r[1,:]))
    a⁻val = zeros(size(r[1,:]))
    κ⁺₁ = zeros(size(r[1,:]))
    κ⁺₂ =zeros(size(r[1,:]))
    κ⁻₁ = zeros(size(r[1,:]))
    κ⁻₂ = zeros(size(r[1,:]))
    F₁⁺₁ = zeros(size(r[1,:]))
    F₁⁺₂ = zeros(size(r[1,:]))
    F₁⁻₁ = zeros(size(r[1,:]))
    F₁⁻₂ =zeros(size(r[1,:]))
    F₂⁺₁ = zeros(size(r[1,:]))
    F₂⁺₂ = zeros(size(r[1,:]))
    F₂⁻₁ = zeros(size(r[1,:]))
    F₂⁻₂ = zeros(size(r[1,:]))
    F₃⁺₁ = zeros(size(r[1,:]))
    F₃⁺₂ = zeros(size(r[1,:]))
    F₃⁻₁ = zeros(size(r[1,:]))
    F₃⁻₂ = zeros(size(r[1,:]))
    H⁺₁ = zeros(size(r[1,:]))
    H⁻₁ = zeros(size(r[1,:]))
    H⁺₂ = zeros(size(r[1,:]))
    H⁻₂ = zeros(size(r[1,:]))
    H⁺₃ = zeros(size(r[1,:]))
    H⁻₃ = zeros(size(r[1,:]))
    L₁ = zeros(size(r[1,:]))
    L₂ = zeros(size(r[1,:]))
    L₃ = zeros(size(r[1,:]))
    ddσ⁺₁ = zeros(size(r[1,:]))
    ddσ⁺₂ = zeros(size(r[1,:]))
    ddσ⁻₁ = zeros(size(r[1,:]))
    ddσ⁻₂ = zeros(size(r[1,:]))

    return r⁺₁,r⁺₂,σ⁺₁,σ⁺₂,η⁺₁,η⁺₂,r⁻₁,r⁻₂,σ⁻₁,σ⁻₂,η⁻₁,η⁻₂,dσ⁺₁,dσ⁺₂,dη⁺₁,dη⁺₂,dσ⁻₁,dσ⁻₂,dη⁻₁,dη⁻₂,β⁺₁,γ⁺₁,β⁺₂,γ⁺₂,λ⁺₁,λ⁺₂ ,a⁺val ,β⁻₁ ,γ⁻₁ ,β⁻₂, γ⁻₂ ,λ⁻₁ ,λ⁻₂, a⁻val ,κ⁺₁ , κ⁺₂ ,κ⁻₁ ,κ⁻₂, F₁⁺₁ ,F₁⁺₂ ,F₁⁻₁ ,F₁⁻₂,F₂⁺₁ ,F₂⁺₂ ,F₂⁻₁,F₂⁻₂,F₃⁺₁,F₃⁺₂,F₃⁻₁, F₃⁻₂, H⁺₁, H⁻₁, H⁺₂, H⁻₂, H⁺₃, H⁻₃, L₁, L₂, L₃, ddσ⁺₁, ddσ⁺₂, ddσ⁻₁, ddσ⁻₂ 
end

function FVM_InitialBoundary(type,r₀,θ,M)
    r = zeros(1,M);

    if type == "circle"
        r .= ones(size(r)).*r₀
    elseif type == "square"
        R = 50 ;
        for i in 1:M
            if θ[i] < π/2
                r[i] = min((R/cos(θ[i])), (R/sin(θ[i])))
            else
                r[i] = min((R/cos((θ[i]%(π/2)))), (R/sin((θ[i]%(π/2)))))
            end
        end
    elseif type =="hex"
        R = √((2/3√3)*π*(r₀^2));
        for i in 1:M
            if θ[i] >= 0 && θ[i] < π/3
                r[i] = R/(cos(θ[i])+1/sqrt(3)*sin(θ[i]))
            elseif θ[i] >= π/3 && θ[i] < π/2
                r[i] = R*sqrt(3)/(2*sin(θ[i]))
            elseif θ[i] >= π/2 && θ[i] < 2*π/3
                r[i] = R*sqrt(3)/(2*sin(θ[i]))
            elseif θ[i] >= 2*π/3 && θ[i] < π
                r[i] = R/(-cos(θ[i])+1/sqrt(3)*sin(θ[i]))
            elseif θ[i] >= π && θ[i] < 4*π/3
                r[i] = -R/(cos(θ[i])+1/sqrt(3)*sin(θ[i]))
            elseif θ[i] >= 4*π/3 && θ[i] < 3*π/2
                r[i] = -R*sqrt(3)/(2*sin(θ[i])) 
            elseif θ[i] >= 3*π/2 && θ[i] < 5*π/3
                r[i] = -R*sqrt(3)/(2*sin(θ[i]))
            else
                r[i] = R/(cos(θ[i])-1/sqrt(3)*sin(θ[i]))
            end
        end
                       
    end
    return r
end