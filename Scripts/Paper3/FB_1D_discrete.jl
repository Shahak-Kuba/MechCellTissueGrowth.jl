import MechCellTissueGrowth as MCTG

Domain = MCTG.DomainProperties_t(N=10, m = 2)
CellMech = MCTG.CellMechProperties_t(kₛ=5, kf = 0)
Prolif = MCTG.CellEvent_t(true, 0.1, "Constant")
Death = MCTG.CellEvent_t()
Embed = MCTG.CellEvent_t()
ProlifEmbed = MCTG.CellEvent_t()
SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)


FB_IC = MCTG.FB_IC_t(q0 = x -> 2.0, q0_der = x -> 0.0, L0 = 10.0)

Seed = 2
NumSaveTimePoints = 1000

sol = MCTG.FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints);

f = Figure(size=(600,600));
ax = Axis(f[1,1], aspect=1, xlabel=L"$x$", ylabel=L"$t$", title=L"$m = 1$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, limits=(-1,50,-1,11));
#lines!(ax, collect(LinRange(0,60, 200)), S.(collect(LinRange(0,60, 200))), color=:black, linewidth=5, label="Substrate");

for jj in 1:200:length(sol.t)
        x = sol.u[jj][1,:]
        t = sol.t[jj]
        #y = sol.u[ii][2,:]
        lines!(ax, x, t*ones(size(x)), color=:black, linewidth=3);
        scatter!(ax, x, t*ones(size(x)), color=:red, markersize=10);
end

display(f)

f = Figure(size=(600,600));
ax = Axis(f[1,1], aspect=1, xlabel=L"$x$", ylabel=L"$q_{i}$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, limits=(0,45,0,3));
#lines!(ax, collect(LinRange(0,60, 200)), S.(collect(LinRange(0,60, 200))), color=:black, linewidth=5, label="Substrate");

for jj in 1:200:length(sol.t)
        x = sol.u[jj][1,:]
        q = sol.Density[jj]
        #y = sol.u[ii][2,:]
        stairs!(ax, x, q, label="t = $(sol.t[jj])", linewidth=3, step=:center);
end

display(f)

"""
m_vals = [2, 4, 10]
all_cell_boundary_positions = []
partial_solutions = []
partial_solutions_t = []
all_solutions_t = []
all_solutions = []
all_solutions_density = []
all_numerical_dl_dt = []
all_final_spring_pos = []

# initial condition (for continuum model)
q0 = x -> 2.0 
q0_der = x -> 0.0
L0 = 10.0

initial_N, error = quadgk(x -> q0(x), 0, L0)

k = 10.0
η = 1.0
a = 2.0

for m in m_vals

    N = round(Int, initial_N); # number of cells
    M = (N * m) + 1; # total number of nodes

    T1 = BBSM.TissueMechProperties_t(k=k, a=a, η=η, restoring_force="hookean")
    TissueMech = BBSM.generate_homogeneous_population(T1, N, m)
    SimTime = BBSM.SimTime_t(Tmax=40.0, δt=0.001, event_δt=0.0)

    # solving problem
    event_cb = PeriodicCallback(TissueRelaxation1D.tissue_dep_affect!, SimTime.event_δt; save_positions = (false, false))
    cbs = CallbackSet(event_cb)

    # substrate function
    S = (x) ->  0 

    # setup initial condition
    y_val = (x) -> 1 #0.1 .* exp.(-((x .- 30).^2) ./ 20)
    u0 = ElasticMatrix{Float64}(undef,2,M)
    u0[1,:] = BBSM.generate_discrete_IC_from_density_profile(q0, q0_der, N, m, L0, 0.01) 
    u0[2,:] = y_val.(u0[1,:])

    prob, p = BBSM.SetupODEproblem(TissueMech, SimTime, u0, S);

    sol = solve(prob, Euler(), dt = SimTime.δt, dtmax = SimTime.δt, callback=cbs);

    all_x_positions = [sol.u[ii][1,1:m:end] for ii in 1:size(sol.u,1)]
    final_spring_pos = [sol.u[ii][1,end-1:end] for ii in 1:size(sol.u,1)]
    density = [1 ./ diff(all_x_positions[ii]) for ii in axes(all_x_positions,1)]

    # density calculation of the boundaries based on Vandenheuvel et al. 2024
    # "Pushing coarse-grain ...."
    for ii in axes(density,1)
        density[ii][1] =  (2 / (all_x_positions[ii][2] - all_x_positions[ii][1])) - (2 / (all_x_positions[ii][3] - all_x_positions[ii][1]))
        density[ii][end] =  (2 / (all_x_positions[ii][end] - all_x_positions[ii][end-1])) - (2 / (all_x_positions[ii][end] - all_x_positions[ii][end-2]))
    end

    boundary_x_positions = [sol.u[ii][1,end] for ii in 1:size(sol.u,1)]
    dl_dt = diff(boundary_x_positions) ./ diff(sol.t[1:size(sol.u,1)])


    push!(all_cell_boundary_positions, all_x_positions)
    push!(partial_solutions, [sol.u[1:10000:end]...])
    push!(partial_solutions_t, sol.t[1:10000:end])
    push!(all_solutions, [sol.u...])
    push!(all_solutions_t, sol.t)
    push!(all_solutions_density, density)
    push!(all_numerical_dl_dt, dl_dt)
    push!(all_final_spring_pos, final_spring_pos)


    println("Completed for m = m")
end
"""