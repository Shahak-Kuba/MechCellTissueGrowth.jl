
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
end

# ODE problem (Spatially discretised)
function rhs!(du, u, p::FBParams, t)
    α, k, a, η, Δx, N = p.α, p.k, p.a, p.η, p.Δx, p.N
    right_BC = p.rBC
    left_BC = p.lBC

    F = p.F
    D = p.D

    diffusivity_method = "arithmetic"

    # Unpack state
    q = u[1:end-1]
    L = u[end]

    q_right_ghost = 0.0

    if right_BC == :fixed
        q_right_ghost = q[N-1]
        dLdt = 0.0
    elseif right_BC == :free
        q_right_ghost = q[N-1] + ((4 * Δx * q[N] * L) / (η * D(q[N]))) * F(q[N]) # ghost node with \mathcal{O}(Δx^2) accuracy
        # to be implemented
        dLdt = (-1/η * F(q[N])) - (D(q[N]) / (2 * q[N] * L)) * ((q_right_ghost - q[N-1]) / (2*Δx))
        #(-2/η) * F(q[N])
        #(-1/η * F(q[N])) - (D(q[N]) / (2 * q[N] * L)) * ((q_right_ghost - q[N-1]) / (2*Δx))
        #(-D(q[N]) / (q[N] * L)) * ((q_right_ghost - q[N-1]) / (2*Δx))
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
            # upwinding on first term (advection term)
            dqidt = (1 / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx) * (Dhp * ( (q_right_ghost - q[ii]) / (Δx) ) - Dhm * ( (q[ii] - q[ii-1]) / (Δx) ) )
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
            dqidt = (z_i / L) * dLdt * ( (q[ii] - q[ii-1])/(Δx) ) + (1/L^2) * (1/Δx) * (Dhp * ( (q[ii+1] - q[ii]) / (Δx) ) - Dhm * ( (q[ii] - q[ii-1]) / (Δx) ) )
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

