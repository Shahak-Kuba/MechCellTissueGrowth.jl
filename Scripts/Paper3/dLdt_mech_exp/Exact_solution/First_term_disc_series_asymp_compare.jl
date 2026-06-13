import MechCellTissueGrowth as MCTG
using CairoMakie
using LambertW
using Printf

# -----------------------------------------------------------------------------
# Parameters (chosen arbitrarily; you can change these)
# -----------------------------------------------------------------------------
N      = 10            # the integer N in theta_p formula
k_star = 5          # k*
eta_star = 1          # eta*  (avoiding name collision with `eta`)
ell0   = 1.0          # l(0)
a_star = 5          # a*
ell0_minus_astar = ell0 - a_star
α  = k_star / eta_star

T = 50

# -----------------------------------------------------------------------------
# Running discrete simulation for m = 5
# -----------------------------------------------------------------------------

m = 20
q₀ = 1.0;
L0 = 10.0;
Domain = MCTG.DomainProperties_t(N=N, m = m, domain_type = "1D")
CellMech = MCTG.CellMechProperties_t(kₛ=k_star, kf = 0, a = a_star, η = eta_star, restoring_force="hookean")
Prolif = MCTG.CellEvent_t()
Death = MCTG.CellEvent_t()
Embed = MCTG.CellEvent_t()
ProlifEmbed = MCTG.CellEvent_t()
SimTime = MCTG.SimTime_t(Tmax=T, δt=0.0001, event_δt=0.0001)
FB_IC = MCTG.FB_IC_t(q0 = x -> q₀, q0_der = x -> 0.0, L0 = L0)
Seed = 1
NumSaveTimePoints = 50001
sol = MCTG.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints);

idx_for_plots = [2, 11, 101, 1001, 10001, 20001, 30001, 40001, 50001]

first_term_disc = -[Domain.m / CellMech.η * CellMech.kₛ * (1/sol.Density[ii][end] - CellMech.a) for ii in eachindex(sol.Density)][idx_for_plots]#[2:4000:end]
t_disc = sol.t[idx_for_plots]#[2:4000:end]

## -----------------------------------------------------------------------------
# Series solution for the first term in dL/dt (continuum limit)
## -----------------------------------------------------------------------------
function L_series(t::Real; P::Int = 200)
    A_inf = (2 / (N * eta_star)) * k_star * ell0_minus_astar
    s = 0.0
    @inbounds for p in 1:P
        s += exp(-α * (2p - 1)^2 * pi^2 * t / (4 * N^2))
    end
    return A_inf * s
end

## -----------------------------------------------------------------------------
# Plotting the series solution and its leading-order short- and long-time asymptotics
## -----------------------------------------------------------------------------

const SERIES_BLACK   = "#000000"   # series form (P=300)
const SHORT_T_BLUE   = "#56B4E9"   # short-t leading order
const LONG_T_RED     = "#D55E00"   # long-t leading order (vermillion)
const POINT_YELLOW   = "#E69F00"   # Lambert-W marker
const STAR_PURPLE    = "#CC79A7"   # discrete m=10 stars

fig2 = Figure(size = (900, 800))
ax2 = Axis(fig2[1, 1];
    xlabel = L"t",
    ylabel = L"|\mathcal{C}_{m}(t)| ",
    xscale = log10,
    yscale = log10,
    xminorticksvisible = true,
    yminorticksvisible = true
)

ts_dense = LinRange(1e-3, T, 400)#10 .^ range(-4, 2; length = 200)
vals_series = [L_series(t; P=300) for t in ts_dense]
#vals_leading_short = [ell0_minus_astar / sqrt(pi * kappa * t) for t in ts_dense]
tau = (k, η, N, t) -> (k*pi^2*t) / (4 * N^2 * η)
vals_leading_short = [(((2 / (N * eta_star)) * k_star * ell0_minus_astar) / 4 ) * sqrt(pi / tau(k_star, eta_star, N, t)) for t in ts_dense]
# Leading-order long-time:   A_inf * exp(-kappa pi^2 t / (4 N^2))
A_inf = (2 / (N * eta_star)) * k_star * ell0_minus_astar
vals_leading_long  = [A_inf * exp(-α * pi^2 * t / (4 * N^2)) for t in ts_dense]
 
lines!(ax2, ts_dense, max.(abs.(vals_series), 1e-16);
    color = SERIES_BLACK, linewidth = 4, label = "series form (P=300)")
#lines!(ax2, ts_dense, max.(abs.(vals_closed), 1e-16);
#    color = :crimson, linewidth = 3, linestyle = :dash, label = "closed form (K=30)")
lines!(ax2, ts_dense, abs.(vals_leading_short);
    color = :blue, linewidth = 4, linestyle = :dash,
    label = "short-t leading order")
lines!(ax2, ts_dense, abs.(vals_leading_long);
    color = :red, linewidth = 4, linestyle = :dash,
    label = "long-t leading order (p=1)")
t_crossover = -lambertw(-π/32, 0) * (8 * N^2 * eta_star)/ (k_star * pi)
scatter!(ax2, t_crossover, abs.(L_series(t_crossover; P=300)); color=:green, markersize=20, label="t ≈ W₀(-π/32)((N² η)/(π k))", strokecolor=:black, strokewidth=0.8)
# plotting discrete solution for m=20
scatter!(ax2, t_disc, abs.(first_term_disc); color=:orange, marker = :star5, markersize=20, label="discrete m = $m", strokecolor=:black, strokewidth=0.8)
axislegend(ax2; position = :lb)
display(fig2)
save("short-long-term-behaviour_p_300_disc_log_log_version.png", fig2)