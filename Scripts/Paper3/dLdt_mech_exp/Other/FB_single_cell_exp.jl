import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK

# single cell simulation to explore δ(t) function effects

m_vals = [5,10,15]
resting_lengths = [2.0,5.0,10.0]
all_solutions = []

for a_value in resting_lengths
    resting_lengths_sols = []
    for m_value in m_vals
        # setting domain size
        L0 = 1.0;
        # solving discrete model
        Domain = MCTG.DomainProperties_t(N=1, m = m_value)
        CellMech = MCTG.CellMechProperties_t(kₛ=10, kf = 0, a = a_value)
        Prolif = MCTG.CellEvent_t()
        Death = MCTG.CellEvent_t()
        Embed = MCTG.CellEvent_t()
        ProlifEmbed = MCTG.CellEvent_t()
        SimTime = MCTG.SimTime_t(Tmax=1, δt=0.000001, event_δt=0.001)

        FB_IC = MCTG.FB_IC_t(q0 = x -> 1.0, q0_der = x -> 0.0, L0 = L0)
        NumSaveTimePoints = 1000

        sol_disc = MCTG.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
        push!(resting_lengths_sols, sol_disc)
    end
    push!(all_solutions, resting_lengths_sols)
end

# plotting evolution of L(t) for different m and a values

f = Figure(size=(1200, 800))
ax1 = Axis(f[1,1], xlabel="Time", ylabel="L(t)", title="a = $(resting_lengths[1])")
ax2 = Axis(f[1,2], xlabel="Time", ylabel="L(t)", title="a = $(resting_lengths[2])")
ax3 = Axis(f[1,3], xlabel="Time", ylabel="L(t)", title="a = $(resting_lengths[3])")
all_axes = [ax1, ax2, ax3]
line_colors = [:blue, :orange, :green]
for i in 1:length(resting_lengths)
    a_value = resting_lengths[i]
    for j in 1:length(m_vals)
        m_value = m_vals[j]
        sol_disc = all_solutions[i][j]
        L_t = [sol_disc.u[ii][1,end] for ii in 1:length(sol_disc.u)]
        t = sol_disc.t
        ax = all_axes[i]
        lc = line_colors[j]
        lines!(ax, t, L_t, label="a=$(a_value), m=$(m_value)", linewidth=3, color=lc)
    end
end
# calculating dL/dt
ax4 = Axis(f[2,1], xlabel="Time", ylabel="L'(t)",  limits=(-0.01,0.05, -20, 1000), xticks=0:0.1:0.25)
ax5 = Axis(f[2,2], xlabel="Time", ylabel="L'(t)",  limits=(-0.01,0.05, -20, 2500), xticks=0:0.1:0.25)
ax6 = Axis(f[2,3], xlabel="Time", ylabel="L'(t)",  limits=(-0.01,0.05, -20, 3000), xticks=0:0.1:0.25)
all_axes = [ax4, ax5, ax6]
line_colors = [:blue, :orange, :green]
for i in 1:length(resting_lengths)
    a_value = resting_lengths[i]
    for j in 1:length(m_vals)
        m_value = m_vals[j]
        sol_disc = all_solutions[i][j]
        L_t = [sol_disc.u[ii][1,end] for ii in 1:length(sol_disc.u)]
        t = sol_disc.t
        dt = t[2] - t[1]
        dL_dt = diff(L_t) ./ dt
        ax = all_axes[i]
        lc = line_colors[j]
        lines!(ax, t[1:end-1], dL_dt, label="m=$(m_value)", linewidth=3, color=lc)
    end
    if i == 3
        axislegend(all_axes[i], position=:rt)
    end
end
display(f)