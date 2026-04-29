# =============================================================================
# verify_limit.jl
#
# Numerically verifies the closed-form limit derived in the LaTeX writeup:
#
#   lim_{m -> inf} m * (l_M(t) - a*) = (l(0) - a*) / sqrt(pi * kappa * t)
#                                    * ( 1 + 2 sum_{j=1}^inf (-1)^j
#                                                  exp(-j^2 N^2 / (kappa t)) )
#
# where kappa = k* / eta*. Equivalently, this is the closed-form value of
# the limiting series A_inf * sum_{p=1}^inf exp(-kappa (2p-1)^2 pi^2 t / (4 N^2)).
#
# Two checks:
#   1. The direct sum (with the corrected A_p) converges to the closed form
#      as m -> inf, at all tested times.
#   2. The series form and the closed (Poisson) form agree wherever both
#      are computable to high precision.
#
# Plots saved with CairoMakie.
#
# Run:  julia --project verify_limit.jl
# =============================================================================
 
using CairoMakie
using Printf
using LambertW
 
# -----------------------------------------------------------------------------
# Parameters (chosen arbitrarily; you can change these)
# -----------------------------------------------------------------------------
N      = 10            # the integer N in theta_p formula
k_star = 5          # k*
eta_st = 1          # eta*  (avoiding name collision with `eta`)
ell0   = 3          # l(0)
a_star = 5          # a*
ell0_minus_astar = ell0 - a_star
kappa  = k_star / eta_st
 
# -----------------------------------------------------------------------------
# Direct evaluation of the original sum (with the corrected A_p).
# This is the LHS of the limit: m * (l_M(t) - a*).
# -----------------------------------------------------------------------------
function L_direct(m::Int, t::Real)
    s = 0.0
    @inbounds for p in 1:(N*m)
        theta = (2p - 1) * pi / (2*N*m + 1)
        lam   = -2 + 2*cos(theta)
        Ap    = (4 * cos(theta/2)^3) / (2*N*m + 1) * m * ell0_minus_astar
        s    += Ap * exp(kappa * m^2 * lam * t)
    end
    return s
end
 
# -----------------------------------------------------------------------------
# Limiting series form (slow at small t).
#   A_inf * sum_{p=1}^P exp(-kappa (2p-1)^2 pi^2 t / (4 N^2))
# -----------------------------------------------------------------------------
function L_series(t::Real; P::Int = 200)
    A_inf = (2 / N) * ell0_minus_astar
    s = 0.0
    @inbounds for p in 1:P
        s += exp(-kappa * (2p - 1)^2 * pi^2 * t / (4 * N^2))
    end
    return A_inf * s
end
 
# -----------------------------------------------------------------------------
# Closed (Poisson) form (fast at small t).
#   (l(0) - a*) / sqrt(pi * kappa * t) *
#       ( 1 + 2 sum_{j=1}^K (-1)^j exp(-j^2 N^2 / (kappa t)) )
# -----------------------------------------------------------------------------
function L_closed(t::Real; K::Int = 20)
    if t <= 0
        return Inf  # singular at t=0
    end
    pref = ell0_minus_astar / sqrt(pi * kappa * t)
    bracket = 1.0
    @inbounds for j in 1:K
        bracket += 2 * (-1)^j * exp(-j^2 * N^2 / (kappa * t))
    end
    return pref * bracket
end
 
# -----------------------------------------------------------------------------
# Adaptive evaluator: chooses series vs. closed form based on t.
# This is the function you would actually use inside the PDE solver.
#
#   tilde_t = kappa * t / N^2  is the dimensionless time.
#   For tilde_t < 1, the closed form converges in O(1) terms (often just K=0).
#   For tilde_t > 1, the series converges in O(1) terms (often just P=1).
# -----------------------------------------------------------------------------
function L_adaptive(t::Real; rtol::Real = 1e-13)
    if t <= 0
        return Inf
    end
    tilde_t = kappa * t / N^2
    if tilde_t < 1.0
        # Use closed form. Estimate K from the smallness of exp(-N^2/(kappa t)).
        # We want exp(-K^2/tilde_t) < rtol  =>  K > sqrt(-tilde_t * log(rtol))
        K = max(1, ceil(Int, sqrt(-tilde_t * log(rtol)) + 2))
        return L_closed(t; K = K)
    else
        # Use series form. We want exp(-(2P-1)^2 pi^2 tilde_t / 4) < rtol.
        # => (2P-1) > (2/pi) sqrt(-log(rtol)/tilde_t)
        P = max(1, ceil(Int, (1/pi) * sqrt(-log(rtol)/tilde_t) + 2))
        return L_series(t; P = P)
    end
end
 
# =============================================================================
# CHECK 1: Convergence of the direct sum to the closed form as m -> inf.
# =============================================================================
 
println("="^78)
println("Check 1: m * (l_M(t) - a*) converges to the closed-form limit as m -> inf.")
println("="^78)
ms = [25, 50, 100, 200, 500, 1000]
ts = [1e-3, 1e-2, 1e-1, 1.0, 5.0]
 
@printf("%10s |", "t")
for m in ms
    @printf("%14s |", "m=$m")
end
@printf("%16s\n", "closed form")
println("-"^(13 + 16*length(ms) + 18))
 
for t in ts
    @printf("%10.3g |", t)
    for m in ms
        @printf("%14.10f |", L_direct(m, t))
    end
    @printf("%16.10f\n", L_closed(t))
end
 
# Save the convergence data for plotting
m_grid    = [10, 20, 50, 100, 200, 500, 1000, 2000]
t_for_conv = [1e-3, 1e-2, 1e-1, 1.0]
errors_vs_m = Dict{Float64, Vector{Float64}}()
for t in t_for_conv
    target = L_closed(t)
    errors_vs_m[t] = [abs(L_direct(m, t) - target) / abs(target) for m in m_grid]
end
 
# =============================================================================
# CHECK 2: Series form vs. Closed form agreement.
# =============================================================================
 
println()
println("="^78)
println("Check 2: Series form (P=300 terms) vs. Closed form (K=30 terms).")
println("="^78)
@printf("%10s | %18s | %18s | %12s\n", "t", "series (P=300)", "closed (K=30)", "rel diff")
println("-"^70)
for t in [1e-6, 1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0]
    s = L_series(t; P=300)
    c = L_closed(t; K=30)
    rd = abs(s-c) / max(abs(s), abs(c), 1e-300)
    @printf("%10.0e | %18.10e | %18.10e | %12.2e\n", t, s, c, rd)
end
 
# =============================================================================
# PLOT 1: Convergence of m * (l_M - a*) to the closed-form limit.
# =============================================================================
 
# Curves: |L_direct(m, t) - L_closed(t)| / |L_closed(t)|  vs. m, for several t.
fig1 = Figure(size = (1200, 1000))
ax1 = Axis(fig1[1, 1];
    xlabel = "m",
    ylabel = "relative error",
    title  = "Convergence of m·(ℓ_M(t) − a*) to the closed-form limit",
    xscale = log10,
    yscale = log10,
)
colors = [:dodgerblue, :crimson, :seagreen, :darkorange]
for (i, t) in enumerate(t_for_conv)
    errs = errors_vs_m[t]
    # protect against floating-point zeros for log scale
    errs_safe = max.(errs, 1e-16)
    scatterlines!(ax1, m_grid, errs_safe; linewidth = 3,
        color = colors[i],
        label = "t = $t",
        markersize = 8,
    )
end
# Reference lines: 1/m and 1/m^2 (the empirical convergence rate)
m_ref = Float64.(m_grid)
lines!(ax1, m_ref, 0.5 ./ m_ref; linewidth = 3,
    color = :gray, linestyle = :dash, label = "1/m (reference)")
axislegend(ax1; position = :rt)
display(fig1)
 
# =============================================================================
# Series vs. Closed form across a wide range of t (log-log).
# =============================================================================
 
fig2 = Figure(size = (900, 800))
ax2 = Axis(fig2[1, 1];
    xlabel = L"t",
    ylabel = L"\log |\mathcal{C}_{m}(t)| ",
    #xscale = log10,
    yscale = log10,
)
 
ts_dense = LinRange(1e-4, 50, 400)#10 .^ range(-4, 2; length = 200)
vals_series = [L_series(t; P=300) for t in ts_dense]
#vals_closed = [L_closed(t; K=30) for t in ts_dense]
# Leading-order short-time:  (l(0)-a*) / sqrt(pi kappa t)
#vals_leading_short = [ell0_minus_astar / sqrt(pi * kappa * t) for t in ts_dense]
tau = (k, η, N, t) -> (k*pi^2*t) / (4 * N^2 * η)
vals_leading_short = [((2 * ell0_minus_astar / N) / 4 ) * sqrt(pi / tau(k_star, eta_st, N, t)) for t in ts_dense]
# Leading-order long-time:   A_inf * exp(-kappa pi^2 t / (4 N^2))
A_inf = (2 / N) * ell0_minus_astar
vals_leading_long  = [A_inf * exp(-kappa * pi^2 * t / (4 * N^2)) for t in ts_dense]
 
lines!(ax2, ts_dense, max.(abs.(vals_series), 1e-16);
    color = :black, linewidth = 4, label = "series form (P=300)")
#lines!(ax2, ts_dense, max.(abs.(vals_closed), 1e-16);
#    color = :crimson, linewidth = 3, linestyle = :dash, label = "closed form (K=30)")
lines!(ax2, ts_dense, abs.(vals_leading_short);
    color = :blue, linewidth = 4, linestyle = :dash,
    label = "short-t leading order")
lines!(ax2, ts_dense, abs.(vals_leading_long);
    color = :red, linewidth = 4, linestyle = :dash,
    label = "long-t leading order (p=1)")
t_crossover = -lambertw(-π/32, 0) * (8 * N^2 * eta_st)/ (k_star * pi)
scatter!(ax2, t_crossover, abs.(L_closed(t_crossover; K=30)); color=:green, markersize=20, label="t ≈ W₀(-π/32)((N² η)/(π k))")
axislegend(ax2; position = :rt)
display(fig2)

save("short-long-term-behaviour_p_300.png", fig2)
 
# =============================================================================
# Number of terms needed in series vs. closed form
# =============================================================================
 
function P_needed_series(t::Real; rtol = 1e-13)
    target = L_series(t; P=500)
    for P in 1:500
        if abs(L_series(t; P=P) - target) / max(abs(target), 1e-300) < rtol
            return P
        end
    end
    return 500
end
 
function K_needed_closed(t::Real; rtol = 1e-13)
    target = L_closed(t; K=50)
    for K in 0:50
        if abs(L_closed(t; K=K) - target) / max(abs(target), 1e-300) < rtol
            return K
        end
    end
    return 50
end
 
ts_for_count = 10 .^ range(-3, 2; length = 60)
P_counts = [P_needed_series(t) for t in ts_for_count]
K_counts = [K_needed_closed(t) for t in ts_for_count]
 
fig3 = Figure(size = (1200, 1000))
ax3 = Axis(fig3[1, 1];
    xlabel = L"\log (t)",
    ylabel = "Terms needed for rel-tol 1e-13",
    title  = "Computational cost: series form vs. closed (Poisson) form",
    xscale = log10,
)
scatterlines!(ax3, ts_for_count, P_counts; linewidth = 3,
    color = :black, markersize = 10, label = "series form (P)")
#scatterlines!(ax3, ts_for_count, K_counts; linewidth = 3,
#    color = :crimson, markersize = 10, label = "closed form (K)")
crossover = N^2 / (kappa * pi)
vlines!(ax3, [crossover];
    color = :gray, linestyle = :dash, linewidth = 3,
    label = "crossover ≈ N²/(κπ) = $(round(crossover, digits=3))")
axislegend(ax3; position = :rt)
display(fig3)