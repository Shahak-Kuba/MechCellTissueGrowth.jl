
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

end

function rhs2!(du, u, p::FBParams, t)
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
