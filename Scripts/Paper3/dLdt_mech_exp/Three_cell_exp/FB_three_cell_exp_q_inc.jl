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
        L0 = 1.1375;
        # solving discrete model
        Domain = MCTG.DomainProperties_t(N=3, m = m_value)
        CellMech = MCTG.CellMechProperties_t(kₛ=5, kf = 0, a = a_value)
        Prolif = MCTG.CellEvent_t()
        Death = MCTG.CellEvent_t()
        Embed = MCTG.CellEvent_t()
        ProlifEmbed = MCTG.CellEvent_t()
        SimTime = MCTG.SimTime_t(Tmax=3, δt=0.0000001, event_δt=0.001)

        #FB_IC = MCTG.FB_IC_t(q0 = x -> 1.0, q0_der = x -> 0.0, L0 = L0)
        FB_IC = MCTG.FB_IC_t(q0 = x -> 2*x .+ 1.5, q0_der = x -> 2.0, L0 = L0)
        NumSaveTimePoints = 1000

        sol_disc = MCTG.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
        push!(resting_lengths_sols, sol_disc)
    end
    push!(all_solutions, resting_lengths_sols)
end

f = Figure(size=(1200, 400))
ax1 = Axis(f[1,1], xlabel="x(t)", ylabel="Time", title="a = $(resting_lengths[1])")
ax2 = Axis(f[1,2], xlabel="x(t)", ylabel="Time", title="a = $(resting_lengths[2])")
ax3 = Axis(f[1,3], xlabel="x(t)", ylabel="Time", title="a = $(resting_lengths[3])")
all_axes = [ax1, ax2, ax3]
for ii in 1:length(resting_lengths)
    sol = all_solutions[ii][1]
    ax = all_axes[ii]
    for jj in 1:50:length(sol.t)
        x = sol.u[jj][1,:]
        t = sol.t[jj] * ones(length(x))
        lines!(ax, x, t, linewidth=3, color=:black)
        scatter!(ax, x, t, markersize=5, color=:red)
    end
end
display(f)



# plotting evolution of L(t) for different m and a values

f2 = Figure(size=(1200, 800))
ax1 = Axis(f2[1,1], xlabel="Time", ylabel="L(t)", title="a = $(resting_lengths[1])")
ax2 = Axis(f2[1,2], xlabel="Time", ylabel="L(t)", title="a = $(resting_lengths[2])")
ax3 = Axis(f2[1,3], xlabel="Time", ylabel="L(t)", title="a = $(resting_lengths[3])")
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
ax4 = Axis(f2[2,1], xlabel="Time", ylabel="L'(t)",  limits=(-0.01,0.05, -20, 1000), xticks=0:0.1:0.25)
ax5 = Axis(f2[2,2], xlabel="Time", ylabel="L'(t)",  limits=(-0.01,0.05, -20, 2500), xticks=0:0.1:0.25)
ax6 = Axis(f2[2,3], xlabel="Time", ylabel="L'(t)",  limits=(-0.01,0.05, -20, 3000), xticks=0:0.1:0.25)
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
display(f2)

# tracking the length of the final cell over time
f3 = Figure(size=(1200, 400))
ax1 = Axis(f3[1,1], xlabel="Time", ylabel="Length", title="a = $(resting_lengths[1])")
ax2 = Axis(f3[1,2], xlabel="Time", ylabel="Length", title="a = $(resting_lengths[2])")
ax3 = Axis(f3[1,3], xlabel="Time", ylabel="Length", title="a = $(resting_lengths[3])")
all_axes = [ax1, ax2, ax3]
line_colors = [:blue, :orange, :green]
for ii in 1:length(resting_lengths)
    for jj in 1:length(m_vals)
        m = m_vals[jj]
        lc = line_colors[jj]
        sol = all_solutions[ii][jj]
        ax = all_axes[ii]
        l_ns = [sol.u[kk][1,end] - sol.u[kk][1,end-m] for kk in 1:length(sol.u)]
        lines!(ax, sol.t, l_ns, linewidth=3,color=lc, label="m=$m")
    end
end
for i in 1:length(resting_lengths)
    a_value = resting_lengths[i]
    for j in 1:length(m_vals)
        m_value = m_vals[j]
        sol_disc = all_solutions[i][j]
        L_t = [sol_disc.u[ii][1,end] for ii in 1:length(sol_disc.u)]
        t = sol_disc.t
        ax = all_axes[i]
        lc = line_colors[j]
        lines!(ax, t, L_t, linestyle=:dash, linewidth=3, color=lc)
    end
end
f3[1, 4] = Legend(f3, ax3, "M values")
display(f3)

f4 = Figure(size=(1200, 800))
ax0 = Axis(f4[1,1:3], xlabel=L"x", title=L"q_{0}(x) = 2x+1.5,\; k/\eta = 5,\; m = 5")
sol = all_solutions[1][1];
x = sol.u[1][1,:];
y = zeros(length(x));
lines!(ax0, x, y, linewidth=5, color=:red)
scatter!(ax0, x, y, markersize=20, color=:black)
ax1 = Axis(f4[2,1], xlabel=L"t", ylabel="Length", title=L"m = 5")
ax2 = Axis(f4[2,2], xlabel=L"t", ylabel="Length", title=L"m = 10")
ax3 = Axis(f4[2,3], xlabel=L"t", ylabel="Length", title=L"m = 15")
all_axes = [ax1, ax2, ax3]
line_colors = [:blue, :orange, :green]
for jj in 1:length(m_vals)
    m = m_vals[jj]
    sol = all_solutions[3][jj]
    N = Int((length(sol.u[1][1,:]) - 1) / m)
    println("N = ", N)
    ax = all_axes[jj]
    cell_lengths = zeros(length(sol.u),N)
    for kk in 1:length(sol.u)
        for nn in 1:N
            cell_lengths[kk,nn] = sol.u[kk][1,nn*m+1] - sol.u[kk][1,(nn-1)*m+1]
        end
    end
    lines!(ax, sol.t, cell_lengths[:,1], linewidth=3,color=line_colors[1], linestyle=:solid, label = "Cell 1")
    lines!(ax, sol.t, cell_lengths[:,2], linewidth=3,color=line_colors[2], linestyle=:solid, label = "Cell 2")
    lines!(ax, sol.t, cell_lengths[:,3], linewidth=3,color=line_colors[3], linestyle=:solid, label = "Cell 3")
    if jj == 3
        axislegend(ax, position=:rb)
    end
end
display(f4)

save("FB_three_cell_Lt_inc.png", f)
save("FB_three_cell_dLdt_inc.png", f2)
save("FB_three_cell_final_cell_length_inc.png", f3)
save("FB_three_cell_all_cell_lengths_inc.png", f4)