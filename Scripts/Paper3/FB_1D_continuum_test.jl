import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK

# setting domain size
L0 = 5.0;

# solving discrete model

Domain = MCTG.DomainProperties_t(N=10, m = 2)
CellMech = MCTG.CellMechProperties_t(kₛ=5, kf = 0)
Prolif = MCTG.CellEvent_t(true, 0.1, "Constant")
Death = MCTG.CellEvent_t()
Embed = MCTG.CellEvent_t()
ProlifEmbed = MCTG.CellEvent_t()
SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)

FB_IC = MCTG.FB_IC_t(q0 = x -> 2.0, q0_der = x -> 0.0, L0 = L0)
NumSaveTimePoints = 1000

# running 100 simulations with different seeds to get an average solution
all_solutions = []
for seed in 1:250
    Seed = seed;
    sol_disc = MCTG.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, seed, NumSaveTimePoints);
    push!(all_solutions, sol_disc)
end

# solution average
avg_solution_cellCount = Float64.(copy(all_solutions[1].CellCount))
for i in 2:length(all_solutions)
    for j in 1:length(avg_solution_cellCount)
        avg_solution_cellCount[j] += all_solutions[i].CellCount[j]
    end
end
avg_solution_cellCount ./= length(all_solutions)


avg_solution_Boundary = zeros(length(all_solutions[1].t))
for i in 1:length(avg_solution_Boundary)
    for j in 1:length(all_solutions)
        avg_solution_Boundary[i] += all_solutions[j].u[i][1,end]
    end
end
avg_solution_Boundary ./= length(all_solutions)


struct FBParams2{TF,TD,TP,TA,TV}
    α::Float64
    k::Float64
    a::Float64
    η::Float64
    Δx::Float64
    L0::Float64
    N::Int
    rBC::Symbol
    lBC::Symbol
    F::TF
    D::TD
    P::TP
    A::TA
    m::Float64
    diffusivity_method::Symbol
    Pq::TV
end

# ------------------------------------------------------------
# Helper: trapezoidal rule on a uniform grid over [0,1]
# vals has length N, spacing Δx = 1/(N-1)
# ------------------------------------------------------------
@inline function trapz_unit_interval(vals, Δx)
    s = 0.5 * vals[1] + 0.5 * vals[end]
    @inbounds for i in 2:length(vals)-1
        s += vals[i]
    end
    return Δx * s
end

# ------------------------------------------------------------
# Optional helper for interface diffusivities
# method = :arithmetic or :harmonic
# ------------------------------------------------------------
@inline function interface_diffusivities(qm, qi, qp, D, method::Symbol)
    Dm = D(qm)
    Di = D(qi)
    Dp = D(qp)

    if method === :arithmetic
        Dhm = 0.5 * (Di + Dm)
        Dhp = 0.5 * (Di + Dp)
    elseif method === :harmonic
        Dhm = (2 * Di * Dm) / (Di + Dm)
        Dhp = (2 * Di * Dp) / (Di + Dp)
    else
        error("Diffusivity method not recognised: $method")
    end

    return Dhm, Dhp
end


function rhs_clean!(du, u, p, t)
    α  = p.α
    η  = p.η
    Δx = p.Δx
    N  = p.N

    right_BC = p.rBC
    left_BC  = p.lBC   # currently unused, kept for symmetry

    F = p.F
    D = p.D
    P = p.P
    A = p.A

    m = p.m
    method = p.diffusivity_method  # :arithmetic or :harmonic

    # State view
    q = @view u[1:N]
    L = u[N+1]

    # Work buffer for P(q)
    Pq = p.Pq

    # Fill Pq in-place to avoid P.(q) allocation
    @inbounds for i in 1:N
        Pq[i] = P(q[i])
    end

    # Cheap quadrature instead of interpolation + quadgk
    g_integral = trapz_unit_interval(Pq, Δx)
    G = L * g_integral

    # Right boundary ghost and L evolution
    q_right_ghost = 0.0
    dLdt = 0.0

    if right_BC === :fixed
        q_right_ghost = q[N-1]
        dLdt = 0.0

    elseif right_BC === :free
        qN   = q[N]
        qNm1 = q[N-1]
        DqN  = D(qN)
        FqN  = F(qN)

        q_right_ghost = qNm1 + ((4 * Δx * qN * L) / DqN) * ((m / η) * FqN + G)

        dLdt = (-(m / η) * FqN) -
               (DqN / (2 * qN * L)) * ((qN - qNm1) / Δx) -
               G
    else
        error("Right BC not recognised: $right_BC")
    end

    invL  = 1 / L
    invL2 = invL * invL
    invΔx = 1 / Δx
    invΔx2 = invΔx * invΔx

    # Left ghost from your original code
    # (Neumann-like symmetry: q_left_ghost = q[2])
    q_left_ghost = q[2]

    # i = 1
    begin
        qi = q[1]
        qp = q[2]
        qm = q_left_ghost

        Dhm, Dhp = interface_diffusivities(qm, qi, qp, D, method)

        reaction = qi * (P(1 / qi) - A(1 / qi))
        diffusion = invL2 * invΔx2 * (Dhp * (qp - qi) - Dhm * (qi - qm))

        du[1] = diffusion + reaction
    end

    # interior nodes: 2:(N-1)
    @inbounds for i in 2:N-1
        qm = q[i-1]
        qi = q[i]
        qp = q[i+1]

        z_i = (i - 1) * Δx

        Dhm, Dhp = interface_diffusivities(qm, qi, qp, D, method)

        # Upwinded advection term, same structure as your original code
        advection = (z_i * invL) * dLdt * ((qi - qm) * invΔx)

        diffusion = invL2 * invΔx *
                    (Dhp * ((qp - qi) * invΔx) -
                     Dhm * ((qi - qm) * invΔx))

        reaction = qi * (P(1 / qi) - A(1 / qi))

        du[i] = advection + diffusion + reaction
    end

    # i = N
    begin
        qm = q[N-1]
        qi = q[N]
        qp = q_right_ghost

        Dhm, Dhp = interface_diffusivities(qm, qi, qp, D, method)

        advection = invL * dLdt * ((qi - qm) * invΔx)

        diffusion = invL2 * invΔx *
                    (Dhp * ((qp - qi) * invΔx) -
                     Dhm * ((qi - qm) * invΔx))

        reaction = qi * (P(1 / qi) - A(1 / qi))

        du[N] = advection + diffusion + reaction
    end

    du[N+1] = dLdt
    return nothing
end

N = 10001
Δx = 1 / (N - 1)

Pq = zeros(N)
# Parameters
k = CellMech.kₛ; η = CellMech.η; α = k/η; a = CellMech.a; N = 10001; m = 1000.0;
# Force and Diffusivity functions (Nonlinear spring model)
Ffunc = ρ -> k * (1/a - ρ)
Dfunc = ρ -> α 
# Proliferation and apoptosis functions
Pfunc = ρ -> 0.1
Afunc = ρ -> 0.0

p = FBParams2(
    α, k, a, η, Δx, L0, N,
    :free,         # rBC
    :fixed,        # lBC
    Ffunc, Dfunc, Pfunc, Afunc,
    m,
    :arithmetic,   # or :harmonic
    Pq
)
q0_func = x -> 2.0 #* cos(z*π/2) + 1.0 #2.0
y0 = MCTG.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)
tspan = (0.0, 10.0)
f = ODEFunction(rhs_clean!)
prob = ODEProblem(f, y0, tspan, p)

sol = solve(
    prob,
    QNDF();
    reltol = 1e-6,
    abstol = 1e-8,
    save_everystep = false
)

using CairoMakie
fig = Figure(size=(650,500))
ax = Axis(fig[1,1], xlabel="x", ylabel="q(x,t)", title="Free boundary")
Nlines = 15
idxs = round.(Int, range(1, length(sol.t), length = min(Nlines, length(sol.t))))
for j in idxs
    x = range(0.0, 1.0, length=p.N) * sol.u[j][end];
    lines!(ax, x, sol.u[j][1:p.N], color=sol.t[j], colorrange = (sol.t[1], sol.t[end]), linewidth=3)
end
Colorbar(fig[1,2], colorrange = (sol.t[1], sol.t[end]), label = "time", width = 15)
display(fig)
#save("qxt_nonlinear_discretisation.png", fig)