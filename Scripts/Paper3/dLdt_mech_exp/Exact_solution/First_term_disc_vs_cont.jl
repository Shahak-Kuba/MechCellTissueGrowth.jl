
import MechCellTissueGrowth as MCTG

m_vals = [1, 5, 10]

f = Figure(size=(600,600))
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"-(m/\eta) * F(q)")

q₀ = 1.0;
N_IC = 10;
k = 5.0; a = 5.0; η = 1.0;

t = LinRange(0, 10, 1001);
F = ρ -> k * (1/ρ - a)
mF_q_term = (2 / (N_IC * η)) .* F(q₀) .* sum(exp.(-(k .* ((2 .* s) - 1).^2 .* π.^2 .* t) ./ (4 .* (N_IC.^2) .* η)) for s in 1:50)
lines!(ax, t, -mF_q_term, linewidth=3, color=:black, label="Asymtotic")


for m in m_vals
    Domain = MCTG.DomainProperties_t(N=N_IC, m = m, domain_type = "1D")
    CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t()
    Embed = MCTG.CellEvent_t()
    ProlifEmbed = MCTG.CellEvent_t()
    SimTime = MCTG.SimTime_t(Tmax=10, δt=0.0001, event_δt=0.0001)
    FB_IC = MCTG.FB_IC_t(q0 = x -> q₀, q0_der = x -> 0.0, L0 = 10.0)
    Seed = 2
    NumSaveTimePoints = 1001
    sol = MCTG.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints);

    # calculating -m/η * F(q) from discrete model
    mF_q_term_disc = -[Domain.m / CellMech.η * CellMech.kₛ * (1/sol.Density[ii][end] - CellMech.a) for ii in eachindex(sol.Density)]
    t = sol.t
    scatter!(ax, 0, mF_q_term_disc[1], markersize=15)
    lines!(ax, t, mF_q_term_disc, linewidth=3, linestyle=:dash, label="Discrete m = $m")
end
axislegend(ax, position=:rt)
display(f)
save("mF_q_comparison.png", f)