function FVM_SolveContinuum_Polar(D,kf,P,A,E,ρ₀,Tmax,r₀,btype, growth_dir)
    # time and space discretisation
    T₀ = 0.0

    N = 10000
    Δt = (Tmax - T₀)/N
    t = LinRange(T₀,Tmax,Int64((Tmax - T₀)/Δt))
 
    m = 241
    Δθ = (2π)/m
    θ = Vector(LinRange(0.0,2π,m))
    pop!(θ)
    M = m - 1

    Φ = 1; # for minmod gradient 

    if growth_dir == "inward"
        S = 1;
    else
        S = -1
    end

    # allocating simulation memory
    ρ = zeros(N+1,M)
    r = zeros(N+1,M)
    dr = zeros(N+1,M)
    ddr = zeros(N+1,M)
    η = zeros(N+1,M)
    dη = zeros(N+1,M)
    σ = zeros(N+1,M)
    dσ = zeros(N+1,M)
    κ = zeros(N+1,M)

    rᵢ = zeros(size(r[1,:]))
    rᵢ₊₁ = zeros(size(r[1,:]))
    rᵢ₋₁ = zeros(size(r[1,:]))
    ηᵢ = zeros(size(η[1,:]))
    ηᵢ₊₁ = zeros(size(η[1,:]))
    ηᵢ₋₁ = zeros(size(η[1,:]))
    σᵢ = zeros(size(σ[1,:]))
    σᵢ₊₁ = zeros(size(σ[1,:]))
    σᵢ₋₁ = zeros(size(σ[1,:]))

    # Preallocating vectors used in simulation (for speed)
    r⁺₁,r⁺₂,σ⁺₁,σ⁺₂,η⁺₁,η⁺₂,r⁻₁,r⁻₂,σ⁻₁,σ⁻₂,
    η⁻₁,η⁻₂,dσ⁺₁,dσ⁺₂,dη⁺₁,dη⁺₂,dσ⁻₁,dσ⁻₂,dη⁻₁,
    dη⁻₂,β⁺₁,Γ⁺₁,β⁺₂,Γ⁺₂,λ⁺₁,λ⁺₂ ,a⁺val ,β⁻₁ ,Γ⁻₁,
    β⁻₂, Γ⁻₂ ,λ⁻₁ ,λ⁻₂, a⁻val ,κ⁺₁ , κ⁺₂ ,κ⁻₁ ,κ⁻₂,
    F₁⁺₁ ,F₁⁺₂ ,F₁⁻₁ ,F₁⁻₂,F₂⁺₁ ,F₂⁺₂ ,F₂⁻₁,F₂⁻₂,F₃⁺₁,
    F₃⁺₂,F₃⁻₁, F₃⁻₂, H⁺₁, H⁻₁, H⁺₂, H⁻₂, H⁺₃, H⁻₃, L₁,L₂,L₃,
    ddσ⁺₁,ddσ⁺₂,ddσ⁻₁,ddσ⁻₂  = SetupSolverMemory(r)


    # setting up initial conditions 
    # (interface: r₀)
    r[1,:] = FVM_InitialBoundary(btype,r₀,θ,M)
    # (σ₀)
    rᵢ .= r[1,:]
    rᵢ₊₁ .= circshift(rᵢ,-1)
    rᵢ₋₁ .= circshift(rᵢ,1)
    σ[1,:] .= (rᵢ₊₁.-rᵢ₋₁)./(2 .*Δθ)

    # (η₀)
    η[1,:] .= ρ₀.*sqrt.(rᵢ.^2 + σ[1,:].^2)
    
    # (ρ₀)
    ρ[1,:] .= ρ₀.*ones(size(ρ[1,:]))


    # simulation time loop 
    @time for n in 1:N

        rᵢ .= r[n,:]
        rᵢ₊₁ .= circshift(rᵢ,-1)
        rᵢ₋₁ .= circshift(rᵢ,1)
        σᵢ .= σ[n,:]
        σᵢ₊₁ .= circshift(σᵢ,-1)
        σᵢ₋₁ .= circshift(σᵢ,1)
        ηᵢ .= η[n,:]
        ηᵢ₊₁ .= circshift(ηᵢ,-1)
        ηᵢ₋₁ .= circshift(ηᵢ,1)

        # getting derivatives
        dr[n,:] .= minmod(rᵢ₋₁, rᵢ, rᵢ₊₁, Φ, Δθ)
        dσ[n,:] .= minmod(σᵢ₋₁, σᵢ, σᵢ₊₁, Φ, Δθ)
        dη[n,:] .= minmod(ηᵢ₋₁, ηᵢ, ηᵢ₊₁, Φ, Δθ)
        ddr[n,:] .= (rᵢ₋₁ .- 2 .*rᵢ .+ rᵢ₊₁)./(Δθ^2)

        # calculate curvature
        κ[n,:] .= -(rᵢ.^2 .+ 2 .* dr[n,:].^2 .- rᵢ.*ddr[n,:])./((rᵢ.^2 .+ dr[n,:].^2).^1.5)
        #κ[n,:] .= -(rᵢ.^2 .- rᵢ.*ddr[n,:] .+ 2 .* dr[n,:])./((rᵢ.*sqrt.(1 .+ (dr[n,:]./rᵢ).^2)).^3)

        # Finite Volume Method Kraganov-Tadmor scheme
        r⁺₁ .= rᵢ₊₁ .- (Δθ/2).*circshift(dr[n,:],-1) 
        r⁺₂ .= rᵢ .+ (Δθ/2).*dr[n,:]
        σ⁺₁ .= σᵢ₊₁ .- (Δθ/2).*circshift(dσ[n,:],-1)
        σ⁺₂ .= σᵢ .+ (Δθ/2).*dσ[n,:]
        η⁺₁ .= ηᵢ₊₁ .- (Δθ/2).*circshift(dη[n,:],-1)
        η⁺₂ .= ηᵢ .+ (Δθ/2).*dη[n,:]

        r⁻₁ .= rᵢ .- (Δθ/2).*dr[n,:]
        r⁻₂ .= rᵢ₋₁ .+ (Δθ/2).*circshift(dr[n,:],1)
        σ⁻₁ .= σᵢ .- (Δθ/2).*dσ[n,:]
        σ⁻₂ .= σᵢ₋₁ .+ (Δθ/2).*circshift(dσ[n,:],1)
        η⁻₁ .= ηᵢ .- (Δθ/2).*dη[n,:]
        η⁻₂ .= ηᵢ₋₁ .+ (Δθ/2).*circshift(dη[n,:],1)

        # derivatives for Kraganov-Tadmor scheme
        dσ⁺₁ .= minmod(circshift(σ⁺₁,1), σ⁺₁, circshift(σ⁺₁,-1), Φ, Δθ)
        dσ⁺₂ .= minmod(circshift(σ⁺₂,1), σ⁺₂, circshift(σ⁺₂,-1), Φ, Δθ)
        dη⁺₁ .= minmod(circshift(η⁺₁,1), η⁺₁, circshift(η⁺₁,-1), Φ, Δθ)
        dη⁺₂ .= minmod(circshift(η⁺₂,1), η⁺₂, circshift(η⁺₂,-1), Φ, Δθ)

        dσ⁻₁ .= minmod(circshift(σ⁻₁,1), σ⁻₁, circshift(σ⁻₁,-1), Φ, Δθ)
        dσ⁻₂ .= minmod(circshift(σ⁻₂,1), σ⁻₂, circshift(σ⁻₂,-1), Φ, Δθ)
        dη⁻₁ .= minmod(circshift(η⁻₁,1), η⁻₁, circshift(η⁻₁,-1), Φ, Δθ)
        dη⁻₂ .= minmod(circshift(η⁻₂,1), η⁻₂, circshift(η⁻₂,-1), Φ, Δθ)

        ddσ⁺₁ .= (circshift(σ⁺₁,1) .- 2 .*σ⁺₁ .+ circshift(σ⁺₁,-1))./(Δθ^2)
        ddσ⁺₂ .= (circshift(σ⁺₂,1) .- 2 .*σ⁺₂ .+ circshift(σ⁺₂,-1))./(Δθ^2)
        ddσ⁻₁ .= (circshift(σ⁻₁,1) .- 2 .*σ⁻₁ .+ circshift(σ⁻₁,-1))./(Δθ^2)
        ddσ⁻₂ .= (circshift(σ⁻₂,1) .- 2 .*σ⁻₂ .+ circshift(σ⁻₂,-1))./(Δθ^2)

        β⁺₁ .= β(r⁺₁, σ⁺₁, η⁺₁, kf, D, dη⁺₁, dσ⁺₁, ddσ⁺₁)
        Γ⁺₁ .= Γ(r⁺₁, σ⁺₁, η⁺₁, kf, D, dσ⁺₁)
        β⁺₂ .= β(r⁺₂, σ⁺₂, η⁺₂, kf, D, dη⁺₂, dσ⁺₂, ddσ⁺₂)
        Γ⁺₂ .= Γ(r⁺₂, σ⁺₂, η⁺₂, kf, D, dσ⁺₂)

        λ⁺₁ .= λ( β⁺₁, Γ⁺₁, r⁺₁, kf)
        λ⁺₂ .= λ( β⁺₂, Γ⁺₂, r⁺₂, kf)
        a⁺val .= a⁺(λ⁺₁, λ⁺₂)

        β⁻₁ .= β(r⁻₁, σ⁻₁, η⁻₁, kf, D, dη⁻₁, dσ⁻₁, ddσ⁻₁)
        Γ⁻₁ .= Γ(r⁻₁, σ⁻₁, η⁻₁, kf, D, dσ⁻₁)
        β⁻₂ .= β(r⁻₂, σ⁻₂, η⁻₂, kf, D, dη⁻₂, dσ⁻₂, ddσ⁻₂)
        Γ⁻₂ .= Γ(r⁻₂, σ⁻₂, η⁻₂, kf, D, dσ⁻₂)

        λ⁻₁ .= λ( β⁻₁, Γ⁻₁, r⁻₁, kf)
        λ⁻₂ .= λ( β⁻₂, Γ⁻₂, r⁻₂, kf)
        a⁻val .= a⁻(λ⁻₁, λ⁻₂)
        
        κ⁺₁ .= -(r⁺₁.^2 .+ 2 .* σ⁺₁.^2 .- r⁺₁.*dσ⁺₁)./((r⁺₁.^2 .+ σ⁺₁.^2).^1.5)
        κ⁺₂ .= -(r⁺₂.^2 .+ 2 .* σ⁺₂.^2 .- r⁺₂.*dσ⁺₂)./((r⁺₂.^2 .+ σ⁺₂.^2).^1.5)
        κ⁻₁ .= -(r⁻₁.^2 .+ 2 .* σ⁻₁.^2 .- r⁻₁.*dσ⁻₁)./((r⁻₁.^2 .+ σ⁻₁.^2).^1.5)
        κ⁻₂ .= -(r⁻₂.^2 .+ 2 .* σ⁻₂.^2 .- r⁻₂.*dσ⁻₂)./((r⁻₂.^2 .+ σ⁻₂.^2).^1.5)

        F₁⁺₁ .= f₁(r⁺₁,σ⁺₁,η⁺₁)
        F₁⁺₂ .= f₁(r⁺₂,σ⁺₂,η⁺₂)
        F₁⁻₁ .= f₁(r⁻₁,σ⁻₁,η⁻₁)
        F₁⁻₂ .= f₁(r⁻₂,σ⁻₂,η⁻₂)

        F₂⁺₁ .= f₂(r⁺₁,σ⁺₁,η⁺₁,kf)
        F₂⁺₂ .= f₂(r⁺₂,σ⁺₂,η⁺₂,kf)
        F₂⁻₁ .= f₂(r⁻₁,σ⁻₁,η⁻₁,kf)
        F₂⁻₂ .= f₂(r⁻₂,σ⁻₂,η⁻₂,kf) 

        F₃⁺₁ .= f₃(r⁺₁,σ⁺₁,η⁺₁,kf,D,dη⁺₁,κ⁺₁)
        F₃⁺₂ .= f₃(r⁺₂,σ⁺₂,η⁺₂,kf,D,dη⁺₂,κ⁺₂)
        F₃⁻₁ .= f₃(r⁻₁,σ⁻₁,η⁻₁,kf,D,dη⁻₁,κ⁻₁)
        F₃⁻₂ .= f₃(r⁻₂,σ⁻₂,η⁻₂,kf,D,dη⁻₂,κ⁻₂)

        H⁺₁ .= 0.5.*(F₁⁺₁ .+ F₁⁺₂) .- ((a⁺val./2).*(r⁺₁ .- r⁺₂))
        H⁻₁ .= 0.5.*(F₁⁻₁ .+ F₁⁻₂) .- ((a⁻val./2).*(r⁻₁ .- r⁻₂))

        H⁺₂ .= 0.5.*(F₂⁺₁ .+ F₂⁺₂) .- ((a⁺val./2).*(σ⁺₁ .- σ⁺₂))
        H⁻₂ .= 0.5.*(F₂⁻₁ .+ F₂⁻₂) .- ((a⁻val./2).*(σ⁻₁ .- σ⁻₂))

        H⁺₃ .= 0.5.*(F₃⁺₁ .+ F₃⁺₂) .- ((a⁺val./2).*(η⁺₁ .- η⁺₂))
        H⁻₃ .= 0.5.*(F₃⁻₁ .+ F₃⁻₂) .- ((a⁻val./2).*(η⁻₁ .- η⁻₂))

        L₁ .= -(1/Δθ).*(H⁺₁ .- H⁻₁) .- S.*kf.*(ηᵢ./rᵢ)
        L₂ .= -(1/Δθ).*(H⁺₂ .- H⁻₂)
        L₃ .= -(1/Δθ).*(H⁺₃ .- H⁻₃) .+ (P - A - E).*ηᵢ

        r[n+1,:] .= r[n,:] + Δt.*L₁
        σ[n+1,:] .= σ[n,:] + Δt.*L₂
        η[n+1,:] .= η[n,:] + Δt.*L₃
        ρ[n+1,:] .= η[n+1,:]./sqrt.(r[n+1,:].^2 + σ[n+1,:].^2)

    end
    return θ,r,ρ,κ,σ,η
end 