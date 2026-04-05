import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK

# setting domain size
L0 = 5.0;

# discrete model Params
CellMech = MCTG.CellMechProperties_t(kₛ=15, kf = 0, a = 1.0, restoring_force="hookean")
Prolif = MCTG.CellEvent_t()
Death = MCTG.CellEvent_t()
Embed = MCTG.CellEvent_t()
ProlifEmbed = MCTG.CellEvent_t()
SimTime = MCTG.SimTime_t(Tmax=100, δt=0.0001, event_δt=0.0001)
FB_IC = MCTG.FB_IC_t(q0 = x -> 1.5, q0_der = x -> 0.0, L0 = L0)
NumSaveTimePoints = 1001

m_vals = [1,2,5,10, 15]

# output from discrete simulations
all_disc_solutions = []
all_disc_solutions_boundary = []
t_disc = []

for m_value in m_vals
    # discrete model Params
    Domain = MCTG.DomainProperties_t(N=10, m = m_value, domain_type="1D")
    # run discrete simulation
    sol_disc = MCTG.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
    push!(all_disc_solutions, sol_disc)
    # data from discrete simulation
    solution_Boundary = zeros(length(sol_disc.t))
    for i in eachindex(sol_disc.t)
        solution_Boundary[i] = sol_disc.u[i][1,end]
    end
    push!(all_disc_solutions_boundary, solution_Boundary)
    push!(t_disc, sol_disc.t)
end

# computing first term
function compute_first_term(disc_sol, m, CellMech)
    t = disc_sol.t
    u = disc_sol.u
    ρ = disc_sol.Density
    a = CellMech.a

    first_term = zeros(size(t))

    for (ii, t_val) in enumerate(t)
        ρ_at_t = ρ[ii]
        first_term[ii] = m * (1 / (m*ρ_at_t[end]) - a/m)
    end
    
    return first_term

end

term_of_interest = []
for (ii, sol) in enumerate(all_disc_solutions)
    first_term = compute_first_term(sol, m_vals[ii], CellMech)
    push!(term_of_interest, first_term)
end

# plotting

f = Figure(size=(600,600));
ax = Axis(f[1,1], aspect=1, xlabel=L"$t$", ylabel=L"$m(\frac{1}{\rho_{M}} - a^{(m)})$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, limits=(-1,20,-0.6,0.1));
for (ii, t_array) in enumerate(t_disc)
    lines!( ax, t_array, term_of_interest[ii], linewidth = 3, label = "m = $(m_vals[ii])" )
end
axislegend(ax, position = :rb)
display(f)

# plot to check power law
f1 = Figure(size=(1000,600));
ax = Axis(f1[1,1], aspect=1, xlabel=L"$m$", ylabel=L"$m(\frac{1}{\rho_{M}} - a^{(m)})$", 
    xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, 
    xscale = log10, yscale = log10);
t_idxs = [1, 101, 201, 301, 501, 601, 1001]
for (ii, t_idx) in enumerate(t_idxs)
    y_val = zeros(size(m_vals))
    for (jj, m_val) in enumerate(m_vals)
        y_val[jj] = abs(term_of_interest[jj][t_idx])
    end
    lines!(ax, m_vals, y_val, linewidth = 3, label = "t = $(t_disc[1][t_idx])" )
    x = log.(m_vals)
    y = log.(y_val)
    ∇ = (y[end] - y[1]) / (x[end] - x[1]) 
    println("Gradient = $∇ at t = $(t_disc[1][t_idx])")
end
Legend(f1[1,2], ax, position = :rb)
display(f1)