import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK
using LinearAlgebra

# setting domain size
L0 = 10.0;
q₀ = 1.5; 

# discrete model Params
CellMech = MCTG.CellMechProperties_t(kₛ=15, kf = 0, a = 1.0, restoring_force="hookean")
Prolif = MCTG.CellEvent_t()
Death = MCTG.CellEvent_t()
Embed = MCTG.CellEvent_t()
ProlifEmbed = MCTG.CellEvent_t()
SimTime = MCTG.SimTime_t(Tmax=50, δt=0.0001, event_δt=0.0001)
FB_IC = MCTG.FB_IC_t(q0 = x -> q₀, q0_der = x -> 0.0, L0 = L0)
NumSaveTimePoints = 1001

m_vals = [1,5,10,20]

# output from discrete simulations
all_disc_solutions = []
all_disc_solutions_boundary = []
t_disc = []

for m_value in m_vals
    # discrete model Params
    Domain = MCTG.DomainProperties_t(N=15, m = m_value, domain_type="1D")
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


# Solving PDE
# Parameters
k = CellMech.kₛ; η = CellMech.η; α = k/η; a = CellMech.a; L0 = L0; N = 1001; Δx = 1/(N-1); 
N_IC = Int(L0 * q₀); rBC=:free; 
tspan = (0.0, SimTime.Tmax);

# Force and Diffusivity functions (Nonlinear spring model)
Ffunc = ρ -> k * (1/ρ - a)
Dfunc = ρ -> α / ρ^2
# Proliferation and apoptosis functions
Pfunc = ρ -> 0.0
Afunc = ρ -> 0.0

q0_func = z -> q₀ 

# solving continuum model without correction term Baker et al. 2019
p = MCTG.FBParams(α=α, η=η, k=k, a=a, L0=L0, N=N, rBC=:free, lBC=:fixed, F=Ffunc, D=Dfunc, P=Pfunc, A=Afunc)
y0 = MCTG.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)
# solve PDE
prob = ODEProblem(MCTG.rhs_Baker!, y0, tspan, p)
sol_Baker = solve(prob, Rodas5P(), saveat=vcat([0.0:0.1:tspan[2]]...))

# solving continuum model with correction term
p = MCTG.FBParams_w_correction(α, k, a, η, Δx, L0, N, N_IC, q₀, rBC, Ffunc, Dfunc, Pfunc, Afunc)
y0 = MCTG.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)
# solve PDE
prob = ODEProblem(MCTG.rhs_with_correction!, y0, tspan, p)
sol = solve(prob, Rodas5P(), saveat=vcat([0.0:0.1:tspan[2]]...))

# plotting L(t) for discrete and continuum models
f = Figure(size=(1300,800));
ax = Axis(f[1,1], aspect=1, xlabel=L"$t$", ylabel=L"$L(t)$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, limits=(-1,51,9,16));
L_cont = [sol_Baker.u[i][end] for i in eachindex(sol_Baker.t)]
lines!(ax, sol_Baker.t, L_cont, label="Baker et al. 2019", linewidth=3);
# continuum solution for L(t) with correction term
L_cont = [sol.u[i][end] for i in eachindex(sol.t)]
lines!(ax, sol.t, L_cont, label="Continuum w/ correction", linewidth=3, color = :black);
# discrete solution for L(t)
clrs = [:red, :green, :orange, :purple]
for ii in eachindex(all_disc_solutions)
    scatter!(ax, all_disc_solutions[ii].t[1:100:end], all_disc_solutions_boundary[ii][1:100:end], label="Discrete m=$(m_vals[ii])", markersize=15, color=clrs[ii]);
    scatter!(ax, all_disc_solutions[ii].t[end], all_disc_solutions_boundary[ii][end], markersize=15, color=clrs[ii]);
end
Legend(f[1,2], ax)
display(f)



"""
f1 = Figure(size=(600,600));
ax = Axis(f1[1,1], aspect=1, xlabel=L"m$", ylabel=L"$||\hat{L}(0) - L^{(m)}(0)||_2$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24)
L2_norm_t_0 = vec(zeros(length(all_disc_solutions),1))
for ii in 1:length(all_disc_solutions)
    L2_norm_t_0[ii] = norm.(L_cont .- all_disc_solutions_boundary[ii])[2]
end
lines!(ax, m_vals, L2_norm_t_0, linewidth=3, color=:black)
display(f1)
"""

f2 = Figure(size=(600,600));
ax = Axis(f2[1,1], aspect=1, xlabel=L"$x$", ylabel=L"$q(x,t)$", title = L"$t = 25$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, limits=(0,15,0.99,1.1));
t_idx = 251;
# continuum solution for L(t)
q_cont = sol.u[t_idx][1:end-1] ;
x = collect(sol.u[t_idx][end] .* range(0, stop=1.0, length=p.N));
lines!(ax, x, q_cont, label="Continuum", linewidth=5);
# discrete solution for L(t)
clrs = [:red, :green, :orange, :purple];
for ii in eachindex(all_disc_solutions)
    dx = diff(all_disc_solutions[ii].u[t_idx][1,:])
    x = (dx[1]/2) * ones(size(all_disc_solutions[ii].u[t_idx][1,1:end-1]));
    x[2:end] = x[2:end] .+ cumsum(dx[2:end])
    stairs!(ax, x, all_disc_solutions[ii].Density[t_idx], label="m=$(m_vals[ii])", linewidth=2, color=clrs[ii]);
end
display(f2)

f3 = Figure(size=(600,600));
ax = Axis(f3[1,1], aspect=1, xlabel=L"$x$", ylabel=L"$q(x,t)$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, limits=(0,10,1.15,1.25));
t_idx = 501;
# continuum solution for L(t)
q_cont = all_cont_solutions[1].u[t_idx][1:end-1] ;
x = collect(all_cont_solutions[1].u[t_idx][end] .* range(0, stop=1.0, length=N));
lines!(ax, x, q_cont, label="Continuum", linewidth=5);
# discrete solution for L(t)
clrs = [:red, :green, :orange, :purple];
for ii in eachindex(all_disc_solutions)
    dx = diff(all_disc_solutions[ii].u[t_idx][1,:])
    x = (dx[1]/2) * ones(size(all_disc_solutions[ii].u[t_idx][1,1:end-1]));
    x[2:end] = x[2:end] .+ cumsum(dx[2:end])
    stairs!(ax, x, all_disc_solutions[ii].Density[t_idx], label="m=$(m_vals[ii])", step=:center, linewidth=2, color=clrs[ii]);
end
#axislegend(ax, position=:rt)

display(f3)


save("disc_vs_cont_fb_Vanden.png", f)
save("disc_vs_cont_den_profile_2.png", f2)
