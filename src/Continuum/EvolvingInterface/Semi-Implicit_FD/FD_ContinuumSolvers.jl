"""
    SolveContinuumLim_Cartesian(D, kf, A, ρ₀, growth_dir, Tmax, Xmax)

Solve the continuum limit of a discrete cell-based model using a semi-implicit scheme and finite differencing for discretization in Cartesian coordinates.

This function is used to simulate the dynamics of a growing cell population in a one-dimensional domain. It utilizes finite difference methods for spatial discretization and a semi-implicit scheme for time integration.

# Arguments
- `D::Float64`: Diffusion coefficient.
- `kf::Float64`: Growth rate constant.
- `A::Float64`: A constant used in the continuum model.
- `ρ₀::Float64`: Initial density of cells.
- `growth_dir::String`: Direction of growth; can be "inward" or "outward".
- `Tmax::Float64`: Maximum time for the simulation.
- `Xmax::Float64`: Maximum spatial extent in the x-direction.

# Returns
- `x::Vector{Float64}`: Vector of spatial positions.
- `h::Matrix{Float64}`: Matrix representing the height profile at each time step and spatial position.
- `ρ::Matrix{Float64}`: Matrix representing the cell density at each time step and spatial position.

# Implementation Details
- The function discretizes the time from `T₀` (initially 0.0) to `Tmax` into `N` steps, and the spatial domain from `X₀` (initially 0.0) to `Xmax` into `M` steps.
- It initializes various matrices and vectors to store derivatives and variables.
- The function sets initial conditions for the height profile `h` and density `ρ`.
- A loop over the time steps computes the finite differences and updates the variables according to the semi-implicit scheme.
- A cyclic tridiagonal matrix method is used to solve a linear system in each time step.
"""
function FD_SolveContinuumLim_Cartesian(D,kf,A,ρ₀,growth_dir,Tmax,Xmax)
    # time and space discretisation
    T₀ = 0.0
    X₀ = 0.0

    N = 3000
    Δt = (Tmax - T₀)/N
    t = LinRange(T₀,Tmax,Int64((Tmax - T₀)/Δt))

    M = 100
    Δx = (Xmax - X₀)/M
    x = LinRange(X₀,Xmax,M)

    # initialise all derivatives and variables
    h = zeros(N+1,M)
    hₓ = zeros(N+1,M)
    hₓₓ = zeros(N+1,M)
    ρ = zeros(N+1,M)
    ρₓ = zeros(N+1,M)
    κ = zeros(N+1,M)
    λ = zeros(N+1,M)
    Φ = zeros(N+1,M)
    Ψ = zeros(N+1,M)

    # initialise matricies for cyclic tridiagonal
    TriDiagMat = zeros(M,M)
    d₀ = zeros(Float64,M)
    d₁ = zeros(Float64,M-1)
    d₋₁ = zeros(Float64,M-1)


    # set initial conditions for h and ρ
    h[1,:] = 2.0 .+ 0.5.*cos.(3*x)
    h[1,end] = h[1,1]
    #x = 3 .* cos.(x)
    ρ[1,:] = ones(1,M)*ρ₀

    if growth_dir == "inward"
        S = -1.0
    else
        S = 1.0
    end

    a⁺ = zeros(size(h[1,:]))
    a⁻ = zeros(size(h[1,:]))
    hᵢ = zeros(size(h[1,:]))
    hᵢ₊₁ = zeros(size(h[1,:]))
    hᵢ₋₁ = zeros(size(h[1,:]))
    ρᵢ = zeros(size(ρ[1,:]))
    ρᵢ₊₁ = zeros(size(ρ[1,:]))
    ρᵢ₋₁ = zeros(size(ρ[1,:]))
    
    @time for n in 1:N

        hᵢ .= h[n,:]
        hᵢ₊₁ .= circshift(hᵢ,-1)
        hᵢ₋₁ .= circshift(hᵢ,1)
        ρᵢ = ρ[n,:]
        ρᵢ₊₁ .= circshift(ρᵢ,-1)
        ρᵢ₋₁ .= circshift(ρᵢ,1)

        a⁺, a⁻ = aₘ(Rᵢ₊₁, Rᵢ₋₁) 

        hₓ[n,:] .= (hᵢ.*(a⁺.-a⁻) .- hᵢ₋₁.*a⁺ .+ hᵢ₊₁.*a⁻) ./ Δx;                    # upwind to find hₓ,  refer eqn (26) in the notes
        hₓₓ[n,:] .= (hᵢ₊₁ .- 2 .*hᵢ .+ hᵢ₋₁) ./ (Δx^2);                               # cental to find hₓₓ, refer eqn (30) in the notes  
        ρₓ[n,:] .= (ρᵢ .* (a⁺.-a⁻) .- ρᵢ₋₁ .* a⁺ .+ ρᵢ₊₁ .* a⁻) ./ Δx;

        κ[n,:] = - hₓₓ[n,:] ./ ( (1 .+ hₓ[n,:].^2).^1.5 )

        Φ[n,:] = ρᵢ .- S.*Δt.*κ[n,:].*kf.*(ρᵢ.^2) .+ (S.*Δt.*kf.*ρᵢ.*ρₓ[n,:].*hₓ[n,:])./(sqrt.(1 .+(hₓ[n,:].^2))) .- D.*Δt.*ρₓ[n,:].*hₓ[n,:].*hₓₓ[n,:] ./ ((1 .+ (hₓ[n,:].^2)).^2) - A.*Δt.*ρᵢ
        λ[n,:] = (D.*Δt./(Δx.^2)) ./ (1 .+ hₓ[n,:].^2)

        h[n+1,:] = hᵢ .- S.*Δt.*kf.*ρᵢ.*sqrt.(1 .+ (hₓ[n,:].^2))

        # cyclic tridiagonal matrix
        d₀ .= 1 .+ 2 .* λ[n,1:M]
        d₁ .= -λ[n,1:M-1]
        d₋₁ .= -λ[n,2:M]
    
        TriDiagMat .= diagm(-1 => d₋₁, 0 => d₀, 1 => d₁)
        TriDiagMat[M,1] = -λ[n,M]
        TriDiagMat[1,M] = -λ[n,1]

        # solving the system
        ρ[n+1,1:M] .= TriDiagMat\Φ[n,1:M]

    end

    return x,h,ρ

end


"""
    SolveContinuumLim_Polar(D, kf, A, ρ₀, growth_dir, Tmax)

Solve the continuum limit of a discrete cell-based model using a semi-implicit scheme and finite differencing for discretization.

This function is specifically designed for a polar coordinate system and is used to model the dynamics of a growing cell population in a circular domain. It employs finite difference methods for spatial discretization and a semi-implicit scheme for time integration.

# Arguments
- `D::Float64`: Diffusion coefficient.
- `kf::Float64`: Growth rate constant.
- `A::Float64`: A constant used in the continuum model.
- `ρ₀::Float64`: Initial density of cells.
- `growth_dir::String`: Direction of growth; can be "inward" or "outward".
- `Tmax::Float64`: Maximum time for the simulation.

# Returns
- `θ::Vector{Float64}`: Vector of angular positions (in radians).
- `R::Matrix{Float64}`: Matrix representing the radial positions at each time step and angular position.
- `ρ::Matrix{Float64}`: Matrix representing the cell density at each time step and angular position.

# Implementation Details
- The function discretizes the time from `T₀` (initially 0.0) to `Tmax` into `N` steps, and the angular domain from 0 to `2π` into `M` steps.
- It initializes various matrices and vectors to store derivatives and variables.
- The function sets initial conditions for radial positions `R` and density `ρ`.
- A loop over the time steps computes the finite differences and updates the variables according to the semi-implicit scheme.
- A cyclic tridiagonal matrix method is used to solve a linear system in each time step.
"""
function FD_SolveContinuumLim_Polar(T,D,kf,A,ρ₀,Tmax,growth_dir,myR)
    # time and space discretisation
    T₀ = 0.0

    N = 10000
    Δt = (Tmax - T₀)/N
    t = LinRange(T₀,Tmax,Int64((Tmax - T₀)/Δt))

    m = 101
    Δθ = (2π)/m
    θ = Vector(LinRange(0.0,2π,m))
    pop!(θ)
    M = m - 1
 

    # initialise all derivatives and variables
    R = zeros(T,N+1,M)
    Rθ = zeros(T,N+1,M)
    Rθθ = zeros(T,N+1,M)
    ρ = zeros(T,N+1,M)
    ρθ = zeros(T,N+1,M)
    κ = zeros(T,N+1,M)
    Φ = zeros(T,N+1,M)
    λ = zeros(T,N+1,M)
    ψ = zeros(T,N+1,M)

    # initialise matricies for cyclic tridiagonalx
    TriDiagMat = zeros(T,M,M)
    d₀ = zeros(T,M)
    d₁ = zeros(T,M-1)
    d₋₁ = zeros(T,M-1)


    # set initial conditions for h and ρ
    #R[1,:] = ones(size(θ)).*3.0 #2.0 .+ 0.5.*cos.(3*θ)
   
    R[1,:] = myR #FD_InitialBoundary(btype,r₀,θ,M)
    ρ[1,:] = ones(T,1,M)*ρ₀

    if growth_dir == "inward"
        S = 1.0
    else
        S = -1.0
    end

    a⁺ = zeros(T,size(R[1,:]))
    a⁻ = zeros(T,size(R[1,:]))
    Rᵢ = zeros(T,size(R[1,:]))
    Rᵢ₊₁ = zeros(T,size(R[1,:]))
    Rᵢ₋₁ = zeros(T,size(R[1,:]))
    ρᵢ = zeros(T,size(ρ[1,:]))
    ρᵢ₊₁ = zeros(T,size(ρ[1,:]))
    ρᵢ₋₁ = zeros(T,size(ρ[1,:]))
    
    for n in 1:N

        Rᵢ .= R[n,:]
        Rᵢ₊₁ .= circshift(Rᵢ,-1)
        Rᵢ₋₁ .= circshift(Rᵢ,1)
        ρᵢ .= ρ[n,:]
        ρᵢ₊₁ .= circshift(ρᵢ,-1)
        ρᵢ₋₁ .= circshift(ρᵢ,1)

        # return Rᵢ₊₁, Rᵢ₋₁
        a⁺, a⁻ = aₘ(Rᵢ₊₁, Rᵢ₋₁) 

        Rθ[n,:] .= (Rᵢ.*(a⁺.-a⁻) .- Rᵢ₋₁.*a⁺ .+ Rᵢ₊₁.*a⁻) ./ Δθ;                        # upwind to find hₓ,  refer eqn (26) in the notes
        Rθθ[n,:] .= (Rᵢ₊₁ .- (2 .*Rᵢ) .+ Rᵢ₋₁) ./ (Δθ^2);                               # cental to find hₓₓ, refer eqn (30) in the notes  
        ρθ[n,:] .= (ρᵢ.*(a⁺.-a⁻) .- ρᵢ₋₁.*a⁺ .+ ρᵢ₊₁.*a⁻) ./ Δθ;

        κ[n,:] .= -(Rᵢ.^2 .+ (2 .* Rθ[n,:].^2) .- (Rᵢ.*Rθθ[n,:]))./((Rᵢ.^2 + Rθ[n,:].^2).^1.5)

        Φ[n,:] .= ρᵢ .- S.*Δt.*κ[n,:].*kf.*(ρᵢ.^2) .- (S.*Δt.*kf.*ρᵢ.*Rθ[n,:].*ρθ[n,:])./(Rᵢ.*sqrt.(Rᵢ.^2 .+ Rθ[n,:].^2)) .- A.*Δt.*ρᵢ
        λ[n,:] .= (D.*Δt./(Δθ.^2)) ./ (Rᵢ.^2 .+ (Rθ[n,:]).^2)
        ψ[n,:] .= (D.*Δt.*Rθ[n,:].*(Rᵢ .+ Rθθ[n,:]))./(Δθ.*((Rᵢ.^2 + Rθ[n,:].^2).^2))

        R[n+1,:] .= Rᵢ .- (S.*Δt.*kf.*ρᵢ./Rᵢ).*sqrt.(Rᵢ.^2 .+ Rθ[n,:].^2)

        # cyclic tridiagonal matrix
        d₀ .= 1 .+ (2 .* λ[n,1:M]) .+ (ψ[n,1:M].*(a⁺ .- a⁻))
        d₁ .= -λ[n,1:M-1] .+ ψ[n,1:M-1].*(a⁻[1:M-1])
        d₋₁ .= -λ[n,2:M] .- ψ[n,2:M].*(a⁺[2:M])
    
        TriDiagMat .= diagm(-1 => d₋₁, 0 => d₀, 1 => d₁)
        TriDiagMat[M,1] = -λ[n,M] + ψ[n,M].*a⁻[M]
        TriDiagMat[1,M] = -λ[n,1] - ψ[n,1].*a⁺[1]

        # solving the system
        ρ[n+1,1:M] .= TriDiagMat\Φ[n,1:M]

    end

    return θ,R,ρ

end

