
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
function rhs!(du, u, p::FBParams, t)
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
        x = range(0, stop=1.0, length=p.N)
        itp = Interpolations.linear_interpolation(x, P.(q));
        g_integral= quadgk(x -> itp(x), 0, 1)[1]
        G = g_integral * L
        #q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (η * D(q[N]))) * F(q[N]) # ghost node with \mathcal{O}(Δx^2) accuracy
        q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (D(q[N]))) * ( (m/η)*F(q[N]) + G )
        dLdt = (-(m/η)* F(q[N])) - (D(q[N])/(2 * q[N] * L))*((q_right_ghost - q[N]) / (Δx))  - G

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

    return nothing
end

function rhs_Baker!(du, u, p::FBParams, t)
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
        m = 100
        x = range(0, stop=1.0, length=p.N)
        itp = Interpolations.linear_interpolation(x, P.(q));
        g_integral= quadgk(x -> itp(x), 0, 1)[1]
        G = g_integral * L
        #q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (η * D(q[N]))) * F(q[N]) # ghost node with \mathcal{O}(Δx^2) accuracy
        q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (D(q[N]))) * ( (1/η)*F(q[N]))
        dLdt = -2/η * F(q[N]) #-(D(q[N])/(2 * q[N] * L))*((q_right_ghost - q[N]) / (Δx))

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
                                   rtol::Real = 1e-12,
                                   t_floor::Real = 1e-10)
    # Add the floor as your code already does, to avoid the 1/√t singularity.
    t_eff = t + t_floor
 
    # Dimensionless time: tilde_t = κ t / N²,  κ = k/η.
    # Crossover at tilde_t ≈ 1/π², but we switch at tilde_t = 1 to be safely
    # within the rapidly-convergent regime of whichever form we pick.
    κ = k / η
    tilde_t = κ * t_eff / N^2
 
    if tilde_t < 1.0
        # ----- Closed (Poisson) form -----
        # mF = (ℓ0 − a) √(k/(π η t)) ( 1 + 2 Σ_{j≥1} (−1)^j exp(−j² N²/(κt)) )
        # Truncation: exp(−K² N²/(κt)) = exp(−K²/tilde_t) < rtol
        # ⇒ K > sqrt(−tilde_t · ln rtol).  Pad by 2 for safety.
        K_max = max(1, ceil(Int, sqrt(-tilde_t * log(rtol)) + 2))
        bracket = 1.0
        @inbounds for j in 1:K_max
            term = 2 * (-1)^j * exp(-j^2 * N^2 / (κ * t_eff))
            bracket += term
            # Early exit if remaining terms are unambiguously below rtol
            if abs(term) < rtol * 1e-2
                break
            end
        end
        return (l0 - a) * sqrt(k / (π * η * t_eff)) * bracket
 
    else
        # ----- Original series form -----
        # mF = (2k/(Nη))(l0 − a) Σ_{p≥1} exp(−k(2p−1)²π² t / (4 N² η))
        #    = A_inf · κ · Σ_{p≥1} exp(−κ(2p−1)²π² t / (4 N²))
        # Truncation: exp(−κ(2P−1)²π² t/(4N²)) = exp(−(2P−1)²π² tilde_t / 4) < rtol
        # ⇒ (2P−1) > (2/π) sqrt(−ln rtol / tilde_t).  Pad by 2 for safety.
        P_max = max(1, ceil(Int, (1/π) * sqrt(-log(rtol) / tilde_t) + 2))
        s = 0.0
        @inbounds for p in 1:P_max
            term = exp(-k * (2p - 1)^2 * π^2 * t_eff / (4 * N^2 * η))
            s += term
            if term < rtol * 1e-2
                break
            end
        end
        return (2 * k / (N * η)) * (l0 - a) * s
    end
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

    #println("t = $t")

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
        mF_q_term = (2 / (N_IC * η)) * (k * (1/q₀ - a)) *
                    sum(exp((-k * (2s - 1)^2 * π^2 * t) / (4 * N_IC^2 * η))
                        for s in 1:50)

        # Linear extrapolation: no artificial flux injected into q dynamics
        q_right_ghost = 2*q[N] - q[N-1] # q[N-1] + ((4 * Δx * q[N] * L) / (D(q[N]))) * ( mF_q_term) # 2*q[N] - q[N-1]

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

    return nothing  # explicit return so the function type is Nothing, not Any
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
