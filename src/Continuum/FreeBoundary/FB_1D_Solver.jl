
@kwdef struct FBParams
    α::Float64 = 1.0
    η::Float64 = 1.0
    k::Float64 = 1.0 
    a::Float64 = 1.0
    N::Int = 101
    L0::Float64 = 10.0
    Δx::Float64 = 1 / (N - 1)
    rBC::Symbol = :fixed  # right boundary condition (:fixed or :free)
    lBC::Symbol = :fixed  # left boundary condition (:fixed or :free)
    F::Function = ρ -> k * (1 / ρ - a) # Force function from Hookean spring
    D::Function = ρ -> α / ρ^2 # Diffusivity function from Hookean spring
    P::Function = ρ -> 0.0 # Proliferation function
    A::Function = ρ -> 0.0 # Apoptosis function
    m::Int = 100 
end

struct FBParams_w_correction{TF, TD, TP, TA}
    α::Float64
    k::Float64
    a::Float64
    η::Float64
    Δx::Float64
    L0::Float64
    N::Int64
    Ncells::Int64
    q₀::Float64
    rBC::Symbol
    F::TF
    D::TD
    P::TP
    A::TA
end

# ODE problem (Spatially discretised)
function rhs_in_ρ!(du, u, p::FBParams, t)
    α, k, a, η, Δx, N = p.α, p.k, p.a, p.η, p.Δx, p.N
    right_BC = p.rBC
    left_BC = p.lBC

    F = p.F
    D = p.D
    P = p.P
    A = p.A


    diffusivity_method = "arithmetic"

    # Unpack state
    q = u[1:end-1]
    L = u[end]

    q_right_ghost = 0.0

    if right_BC == :fixed
        q_right_ghost = q[N-1]
        dLdt = 0.0
    elseif right_BC == :free
        m = 1
        #x = range(0, stop=1.0, length=p.N)
        #itp = Interpolations.linear_interpolation(x, P.(q));
        #g_integral= quadgk(x -> itp(x), 0, 1)[1]
        #G = g_integral * L
        #q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (η * D(q[N]))) * F(q[N]) # ghost node with \mathcal{O}(Δx^2) accuracy
        q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (D(q[N]))) * ( (m/η)*F(q[N]) )#+ G )
        dLdt = (-(1/η)* F(q[N])) - (D(q[N])/(2 * q[N] * L))*((q_right_ghost - q[N]) / (Δx))  #- G

    end

    dqidt = 0.0
    for ii in 1:N
        if ii == 1
            q_left_ghost = q[2]
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q_left_ghost) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q_left_ghost)) / (D(q[ii]) + D(q_left_ghost))
            else
                error("Diffusivity method not recognised")
            end
            dqidt = (1/L^2) * (1/Δx^2) * (Dhp * ( (q[ii+1] - q[ii]) ) - Dhm * ( (q[ii] - q_left_ghost) ) ) + q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        elseif ii == N
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q_right_ghost)
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q_right_ghost)) / (D(q[ii]) + D(q_right_ghost))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (1 / L)*dLdt*( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * ( Dhp*(q_right_ghost - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) + q[ii]*( P(1/q[ii]) - A(1/q[ii]) )
        else
            z_i = (ii-1) * Δx
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (z_i / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * (Dhp*(q[ii+1] - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) + q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        end
        du[ii] = dqidt
    end

    du[N + 1] = dLdt

    return nothing
end

function rhs_Baker_in_ρ!(du, u, p::FBParams, t)
    α, k, a, η, Δx, N = p.α, p.k, p.a, p.η, p.Δx, p.N
    right_BC = p.rBC
    left_BC = p.lBC

    F = p.F
    D = p.D
    P = p.P
    A = p.A

    m = p.m

    diffusivity_method = "arithmetic"

    # Unpack state
    q = u[1:end-1]
    L = u[end]

    q_right_ghost = 0.0

    if right_BC == :fixed
        q_right_ghost = q[N-1]
        dLdt = 0.0
    elseif right_BC == :free
        #m = 100
        #x = range(0, stop=1.0, length=p.N)
        #itp = Interpolations.linear_interpolation(x, P.(q));
        #g_integral= quadgk(x -> itp(x), 0, 1)[1]
        #G = g_integral * L
        #q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (η * D(q[N]))) * F(q[N]) # ghost node with \mathcal{O}(Δx^2) accuracy
        q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (D(q[N]))) * ( (1/η)*F(q[N]))
        dLdt = -(D(q[N])/(q[N] * L))*((q_right_ghost - q[N-1]) / (2*Δx)) #-2/η * F(q[N]) #-(D(q[N])/(2 * q[N] * L))*((q_right_ghost - q[N]) / (Δx))

    end

    dqidt = 0.0
    for ii in 1:N
        if ii == 1
            #dqidt = ( (2 * α)/(L^2 * q[1]^2) ) * ( ((q[2] - q[1]) / p.Δx^2) )
            q_left_ghost = q[2]
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q_left_ghost) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q_left_ghost)) / (D(q[ii]) + D(q_left_ghost))
            else
                error("Diffusivity method not recognised")
            end
            dqidt = (1/L^2) * (1/Δx^2) * (Dhp * ( (q[ii+1] - q[ii]) ) - Dhm * ( (q[ii] - q_left_ghost) ) ) + q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        elseif ii == N
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q_right_ghost)
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q_right_ghost)) / (D(q[ii]) + D(q_right_ghost))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (1 / L)*dLdt*( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * ( Dhp*(q_right_ghost - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) + q[ii]*( P(1/q[ii]) - A(1/q[ii]) )
        else
            z_i = (ii-1) * Δx
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (z_i / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * (Dhp*(q[ii+1] - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) + q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        end
        du[ii] = dqidt
    end

    du[N + 1] = dLdt

end

@inline function compute_mF_q_term(t::Real, k::Real, η::Real, N::Integer,
                                   l0::Real, a::Real;
                                   t_floor::Real = 1e-10,
                                   τ::Real = 1.1303809932021487)
    # Add the floor to avoid the 1/√t singularity.
    t_eff = t + t_floor
    
    t_switch = (2 * N^2 * η) / (k * π^2) * τ

    A_tilde = (2/N) * (l0 - a)
    τ_at_t = (k*π^2*t_eff)/(4 * N^2 * η)
    
    if t_eff < t_switch
        return (A_tilde / 2) * √(pi/τ_at_t)
    else
        return A_tilde * exp(-(k * 0.25 * pi^2 * t_eff)/(2 * N^2 * η))
    end
end

function L_series(t::Real, k::Real, η::Real, N::Integer, l0::Real, a::Real; P::Int = 200)
    A_inf = (2 / (N * η)) * k * (l0 - a)
    α = k / η
    s = 0.0
    @inbounds for p in 1:P
        s += exp(-α * (2p - 1)^2 * pi^2 * t / (4 * N^2))
    end
    return A_inf * s
end


function rhs_Baker_in_q!(du, u, p::FBParams_w_correction, t)
    α, k, a, η, Δx, N = p.α, p.k, p.a, p.η, p.Δx, p.N
    N_IC = p.Ncells
    q₀    = p.q₀
    right_BC = p.rBC

    F = p.F
    D = p.D
    P = p.P
    A = p.A

    diffusivity_method = "arithmetic"

    # Unpack state
    q = u[1:end-1]
    L = u[end]

    q_right_ghost = 0.0

    if right_BC == :fixed
        q_right_ghost = q[N-1]
        dLdt = 0.0
    elseif right_BC == :free
        mF_q_term = L_series(t,k, η, N_IC, 1/q₀, a) #compute_mF_q_term(t,k, η, N_IC, 1/q₀, a)
                    #(2 / (N_IC * η)) * (k * (1/q₀ - a)) *
                    #sum(exp((-k * (2s - 1)^2 * π^2 * t) / (4 * N_IC^2 * η))
                    #    for s in 1:50)

        q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (D(q[N]))) * ( mF_q_term )
        dLdt = -2 * mF_q_term #-(D(q[N])/(q[N] * L))*((q_right_ghost - q[N-1]) / (2Δx)) #-2/η * F(q[N]) #-(D(q[N])/(2 * q[N] * L))*((q_right_ghost - q[N]) / (Δx))

    end

    dqidt = 0.0
    for ii in 1:N
        if ii == 1
            #dqidt = ( (2 * α)/(L^2 * q[1]^2) ) * ( ((q[2] - q[1]) / p.Δx^2) )
            q_left_ghost = q[2]
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q_left_ghost) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q_left_ghost)) / (D(q[ii]) + D(q_left_ghost))
            else
                error("Diffusivity method not recognised")
            end
            dqidt = (1/L^2) * (1/Δx^2) * (Dhp * ( (q[ii+1] - q[ii]) ) - Dhm * ( (q[ii] - q_left_ghost) ) ) #+ q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        elseif ii == N
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q_right_ghost)
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q_right_ghost)) / (D(q[ii]) + D(q_right_ghost))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (1 / L)*dLdt*( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * ( Dhp*(q_right_ghost - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) #+ q[ii]*( P(1/q[ii]) - A(1/q[ii]) )
        else
            z_i = (ii-1) * Δx
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (z_i / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * (Dhp*(q[ii+1] - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) #+ q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        end
        du[ii] = dqidt
    end

    du[N + 1] = dLdt

end

function rhs_with_correction!(du, u, p::FBParams_w_correction, t)
    α, k, a, η, Δx, N = p.α, p.k, p.a, p.η, p.Δx, p.N
    N_IC = p.Ncells
    q₀    = p.q₀
    right_BC = p.rBC

    F = p.F
    D = p.D
    P = p.P
    A = p.A

    diffusivity_method = "arithmetic"

    # Unpack state
    q = u[1:end-1]
    L = u[end]

    # Initialise before conditional so Julia can infer concrete types
    q_right_ghost = 0.0
    dLdt          = 0.0

    if right_BC == :fixed
        q_right_ghost = q[N-1]
        dLdt          = 0.0

    elseif right_BC == :free
        #mF_q_term = compute_mF_q_term(t, k, η, N, 1/q₀, a; rtol=1e-12, t_floor=1e-10)
        mF_q_term = L_series(t, k, η, N_IC, 1/q₀, a) #compute_mF_q_term(t,k, η, N_IC, 1/q₀, a)
                    #(2 / (N_IC * η)) * (k * (1/q₀ - a)) *
                    #sum(exp((-k * (2s - 1)^2 * π^2 * t) / (4 * N_IC^2 * η))
                    #    for s in 1:50)

        # Linear extrapolation: no artificial flux injected into q dynamics
        q_right_ghost = 2*q[N] - q[N-1]
        #q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (D(q[N]))) * ( mF_q_term )

        # dLdt from the asymptotic formula + actual gradient
        # (For uniform IC the gradient term is ~0; mF_q_term carries all the physics)
        dLdt = -mF_q_term - (D(q[N]) / (2 * q[N] * L)) * (q[N] - q[N-1]) / Δx
    end

    dqidt = 0.0

    for ii in 1:N
        if ii == 1
            #dqidt = ( (2 * α)/(L^2 * q[1]^2) ) * ( ((q[2] - q[1]) / p.Δx^2) )
            q_left_ghost = q[2]
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q_left_ghost) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q_left_ghost)) / (D(q[ii]) + D(q_left_ghost))
            else
                error("Diffusivity method not recognised")
            end
            dqidt = (1/L^2) * (1/Δx^2) * (Dhp * ( (q[ii+1] - q[ii]) ) - Dhm * ( (q[ii] - q_left_ghost) ) ) #+ q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        elseif ii == N
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q_right_ghost)
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q_right_ghost)) / (D(q[ii]) + D(q_right_ghost))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (1 / L)*dLdt*( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * ( Dhp*(q_right_ghost - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) #+ q[ii]*( P(1/q[ii]) - A(1/q[ii]) )
        else
            z_i = (ii-1) * Δx
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (z_i / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx^2) * (Dhp*(q[ii+1] - q[ii]) - Dhm*(q[ii] - q[ii-1]) ) #+ q[ii] * ( P(1/q[ii]) - A(1/q[ii]) )
        end
        du[ii] = dqidt
    end

    du[N + 1] = dLdt

    return nothing  # explicit return so the function type is Nothing, not Any
end

function rhs_with_correction!(du, u, p::FBParams_w_correction, t)
    α, k, a, η, Δx, N = p.α, p.k, p.a, p.η, p.Δx, p.N
    N_IC = p.Ncells
    q₀    = p.q₀
    right_BC = p.rBC

    F = p.F
    D = p.D
    P = p.P
    A = p.A

    diffusivity_method = "arithmetic"

    # Unpack state
    q = u[1:end-1]
    L = u[end]

    # Initialise before conditional so Julia can infer concrete types
    q_right_ghost = 0.0
    dLdt          = 0.0

    if right_BC == :fixed
        q_right_ghost = q[N-1]
        dLdt          = 0.0

    elseif right_BC == :free
        # m/η* · F̃(q) in the m→∞ limit (analytical, always ≥ 0 here)
        mF_q_term = compute_mF_q_term(t, k, η, N_IC, 1/q₀, a)

        # --- CHANGE 1: physical Robin ghost node -------------------------------
        # BC at x=L:  -mF̃ + D(q)/(2q) ∂q/∂x = 0  ⟹  D(q)/(2q) ∂q/∂x = mF_q_term
        # On the ξ-grid (x = Lξ):  ∂q/∂x = (1/L)(q_ghost - q[N-1])/(2Δx)
        # Solving for the ghost:
        #   q_ghost = q[N-1] + (4 Δx L q[N] / D(q[N])) · mF_q_term
        # This is what actually feeds the boundary flux into q. The old
        # linear extrapolation (2q[N]-q[N-1]) collapsed the boundary flux to
        # the interior value and silently imposed a no-flux-divergence BC,
        # which is why q stayed frozen.
        q_right_ghost = q[N-1] + (4 * Δx * L * q[N] / D(q[N])) * mF_q_term

        # --- CHANGE 2: BC-consistent dLdt --------------------------------------
        # dL/dt = [-mF̃ - D(q)/(2q) ∂q/∂x] at x=L.
        # Use the SAME centred gradient as the ghost (q_ghost - q[N-1])/(2Δx),
        # and the SAME conversion ∂q/∂x = (1/L)·∂q/∂ξ. By construction this
        # reduces to dLdt = -2·mF_q_term, but we write it explicitly so it
        # stays correct if the BC/ghost is ever changed.
        dqdx_L = (1 / L) * (q_right_ghost - q[N-1]) / (2 * Δx)
        dLdt   = -mF_q_term - (D(q[N]) / (2 * q[N])) * dqdx_L
    end

    dqidt = 0.0

    for ii in 1:N
        if ii == 1
            q_left_ghost = q[2]
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q_left_ghost) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q_left_ghost)) / (D(q[ii]) + D(q_left_ghost))
            else
                error("Diffusivity method not recognised")
            end
            dqidt = (1/L^2) * (1/Δx^2) * (Dhp * ( (q[ii+1] - q[ii]) ) - Dhm * ( (q[ii] - q_left_ghost) ) )
        elseif ii == N
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q_right_ghost)
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q_right_ghost)) / (D(q[ii]) + D(q_right_ghost))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # --- CHANGE 3: mesh-motion term at the boundary node ---------------
            # z_N = (N-1)*Δx must equal 1 on a [0,1] ξ-grid. Use z_N explicitly
            # so the coefficient is consistent with the interior nodes
            # regardless of how Δx is defined.
            z_N = (N - 1) * Δx
            # upwinding on first term (advection term)
            dqidt = (z_N / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) +
                    (1/L^2) * (1/Δx^2) * ( Dhp*(q_right_ghost - q[ii]) - Dhm*(q[ii] - q[ii-1]) )
        else
            z_i = (ii-1) * Δx
            # different diffusivity averaging methods
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1]) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (advection term)
            dqidt = (z_i / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) +
                    (1/L^2) * (1/Δx^2) * (Dhp*(q[ii+1] - q[ii]) - Dhm*(q[ii] - q[ii-1]) )
        end
        du[ii] = dqidt
    end

    du[N + 1] = dLdt

    return nothing  # explicit return so the function type is Nothing, not Any
end

function rhs_with_correction!(du, u, p::FBParams_w_correction, t)
    α, k, a, η, Δx, N = p.α, p.k, p.a, p.η, p.Δx, p.N
    N_IC = p.Ncells
    q₀    = p.q₀
    right_BC = p.rBC

    F = p.F
    D = p.D
    P = p.P
    A = p.A

    diffusivity_method = "arithmetic"

    # Unpack state
    q = u[1:end-1]
    L = u[end]

    # Initialise before conditional so Julia can infer concrete types
    q_right_ghost = 0.0
    dLdt          = 0.0

    if right_BC == :fixed
        q_right_ghost = q[N-1]
        dLdt          = 0.0

    elseif right_BC == :free
        mF_q_term = L_series(t, k, η, N_IC, 1/q₀, a)
        dqdξ_N    = (q[N] - q[N-1]) / Δx          # one-sided, on the ξ-grid
        dLdt      = -mF_q_term - (D(q[N]) / (2 * q[N] * L)) * dqdξ_N
    end

    dqidt = 0.0

    for ii in 1:N
        if ii == 1
            q_left_ghost = q[2]
            if diffusivity_method == "arithmetic"
                Dm = D(q_left_ghost) # central finite difference ghost node
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q_left_ghost)) / (D(q[ii]) + D(q_left_ghost))
            else
                error("Diffusivity method not recognised")
            end
            dqidt = (1/L^2) * (1/Δx^2) * (Dhp * ( (q[ii+1] - q[ii]) ) - Dhm * ( (q[ii] - q_left_ghost) ) )
            du[ii] = dqidt

        elseif ii == N
            # m → ∞: Dirichlet boundary, q[N] is pinned to 1/a.
            # Do NOT integrate this node — freeze it so the constraint holds.
            if right_BC == :free
                du[ii] = 0.0
            else
                # :fixed branch keeps the original interior-style update
                if diffusivity_method == "arithmetic"
                    Dm = D(q[ii-1])
                    Di = D(q[ii])
                    Dp = D(q_right_ghost)
                    Dhp = 0.5 * (Di + Dp)
                    Dhm = 0.5 * (Di + Dm)
                elseif diffusivity_method == "harmonic"
                    Dhp = (2*D(q[ii])*D(q_right_ghost)) / (D(q[ii]) + D(q_right_ghost))
                    Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
                else
                    error("Diffusivity method not recognised")
                end
                z_N = (N - 1) * Δx
                dqidt = (z_N / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) +
                        (1/L^2) * (1/Δx^2) * ( Dhp*(q_right_ghost - q[ii]) - Dhm*(q[ii] - q[ii-1]) )
                du[ii] = 0#dqidt
            end

        else
            z_i = (ii-1) * Δx
            if diffusivity_method == "arithmetic"
                Dm = D(q[ii-1])
                Di = D(q[ii])
                Dp = D(q[ii+1])
                Dhp = 0.5 * (Di + Dp)
                Dhm = 0.5 * (Di + Dm)
            elseif diffusivity_method == "harmonic"
                Dhp = (2*D(q[ii])*D(q[ii+1])) / (D(q[ii]) + D(q[ii+1]))
                Dhm = (2*D(q[ii])*D(q[ii-1])) / (D(q[ii]) + D(q[ii-1]))
            else
                error("Diffusivity method not recognised")
            end
            # upwinding on first term (mesh-motion / advection term)
            dqidt = (z_i / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) +
                    (1/L^2) * (1/Δx^2) * (Dhp*(q[ii+1] - q[ii]) - Dhm*(q[ii] - q[ii-1]) )
            du[ii] = dqidt
        end
    end

    du[N + 1] = dLdt

    return nothing
end

# Build an initial condition U0(z) on the grid
function make_initial_condition_FB(N; U0fun = z -> 1.0 , L0 = 5.0)
    z = range(0.0, 1.0, length=N)
    U0 = [U0fun(zi) for zi in z]
    return vcat(U0, L0)
end

function FBPDEsolve(p::FBParams, y0, tspan) 
    prob = ODEProblem(rhs!, y0, tspan, p)
    sol = solve(prob, Rodas5P(), saveat=vcat([0.0:0.0001:tspan[2]]...))
    return sol
end
