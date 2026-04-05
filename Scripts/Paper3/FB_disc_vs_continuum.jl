import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK

# setting domain size
L0 = 10.0;

# solving discrete model

Domain = MCTG.DomainProperties_t(N=15, m = 5, domain_type = "1D")
CellMech = MCTG.CellMechProperties_t(kₛ=6.7, kf = 0, a = 1.0, restoring_force="hookean")
Prolif = MCTG.CellEvent_t(true, 0.1, "Constant")
Death = MCTG.CellEvent_t()
Embed = MCTG.CellEvent_t()
ProlifEmbed = MCTG.CellEvent_t()
SimTime = MCTG.SimTime_t(Tmax=30, δt=0.001, event_δt=0.001)

FB_IC = MCTG.FB_IC_t(q0 = x -> 1.5, q0_der = x -> 0.0, L0 = L0)
NumSaveTimePoints = 1001

# running 100 simulations with different seeds to get an average solution
all_solutions = []
for seed in 1:500
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

# computing errors
upper_error_cellCount = zeros(size(avg_solution_cellCount))
lower_error_cellCount = zeros(size(avg_solution_cellCount))

for tt in eachindex(all_solutions[1].t) # looping over time
    sum_of_sq_diff = 0.0
    for ii in eachindex(all_solutions)
        sum_of_sq_diff = sum_of_sq_diff + (avg_solution_cellCount[tt] - all_solutions[ii].CellCount[tt])^2
    end
    println(sum_of_sq_diff)
    σ = sqrt((1/(length(all_solutions) - 1)) * sum_of_sq_diff)
    SE = σ / sqrt(length(all_solutions))
    # at a 95% confidence interval
    upper_error_cellCount[tt] = SE
    lower_error_cellCount[tt] = SE
end


avg_solution_Boundary = zeros(length(all_solutions[1].t))
for i in 1:length(avg_solution_Boundary)
    for j in 1:length(all_solutions)
        avg_solution_Boundary[i] += all_solutions[j].u[i][1,end]
    end
end
avg_solution_Boundary ./= length(all_solutions)

# computing errors
upper_error_Boundary= zeros(size(avg_solution_Boundary))
lower_error_Boundary = zeros(size(avg_solution_Boundary))

for tt in eachindex(all_solutions[1].t) # looping over time
    sum_of_sq_diff = 0.0
    for ii in eachindex(all_solutions)
        sum_of_sq_diff = sum_of_sq_diff + (avg_solution_Boundary[tt] - all_solutions[ii].u[tt][1,end])^2
    end
    #println(sum_of_sq_diff)
    σ = sqrt((1/(length(all_solutions) - 1)) * sum_of_sq_diff)
    SE = σ / sqrt(length(all_solutions))
    # at a 95% confidence interval
    upper_error_Boundary[tt] = SE
    lower_error_Boundary[tt] = SE
end

# solving continuum model
# Parameters
k = CellMech.kₛ; η = CellMech.η; α = k/η; a = CellMech.a; L0 = L0; N = 1001;
# Force and Diffusivity functions (Nonlinear spring model)
Ffunc = ρ -> k * (1/ρ - a)
Dfunc = ρ -> α / ρ^2
# Proliferation and apoptosis functions
Pfunc = ρ -> 0.1
Afunc = ρ -> 0.0

p = MCTG.FBParams(α=α, η=η, k=k, a=a, L0=L0, N=N, rBC=:free, lBC=:fixed, F=Ffunc, D=Dfunc, P=Pfunc, A=Afunc, m = 10)
q0_func = x -> 1.5 #* cos(z*π/2) + 1.0 #2.0
y0 = MCTG.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)
tspan = (0.0, SimTime.Tmax)
# solve PDE
prob = ODEProblem(MCTG.rhs!, y0, tspan, p)
#sol_cont = solve(prob, QNDF(); reltol=1e-6, abstol=1e-8, saveat=vcat([0.0:0.01:SimTime.Tmax]...));
sol_cont = solve(prob, Rodas5P(), saveat=vcat([0.0:0.01:SimTime.Tmax]...));


# plotting L(t) for discrete and continuum models
f = Figure(size=(600,600));
ax = Axis(f[1,1], aspect=1, xlabel=L"$t$", ylabel=L"$L(t)$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, limits=(-2,32,10,30));
# continuum solution for L(t)
L_cont = [sol_cont.u[i][end] for i in eachindex(sol_cont.t)]
lines!(ax, sol_cont.t, L_cont, label="Continuum", linewidth=3, color=:blue);
# discrete solution for L(t)
scatter!(ax, all_solutions[1].t[1:100:end], avg_solution_Boundary[1:100:end], label="Discrete", markersize=10, color=:red);
scatter!(ax, all_solutions[1].t[end], avg_solution_Boundary[end], label="Discrete", markersize=10, color=:red);
errorbars!(ax, all_solutions[1].t[1:100:end], avg_solution_Boundary[1:100:end], lower_error_Boundary[1:100:end], upper_error_Boundary[1:100:end], whiskerwidth = 10, color = :red)
display(f)


# plotting Number of cells vs time for discrete and continuum models
f2 = Figure(size=(600,600));
ax = Axis(f2[1,1], aspect=1, xlabel=L"$t$", ylabel=L"$N(t)$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, limits=(0,31,-10,310));
# discrete solution for cell count
t = all_solutions[1].t[1:50:end]
N = avg_solution_cellCount[1:50:end]
# continuum solution for cell count
using Interpolations
cont_cell_count = []
for ii in eachindex(sol_cont.t)
    q = sol_cont.u[ii][1:end-1]
    x = sol_cont.u[ii][end] .* range(0, stop=1.0, length=p.N)
    itp = Interpolations.linear_interpolation(x, q);
    integral, error = quadgk(x -> itp(x), 0, sol_cont.u[ii][end])
    push!(cont_cell_count, integral)
end
# plotting
lines!(ax, sol_cont.t, cont_cell_count, label="Continuum", linewidth=3, color=:blue);
errorbars!(ax, t, N, lower_error_cellCount[1:50:end], upper_error_cellCount[1:50:end], whiskerwidth = 10, color = :red)
scatter!(ax, t, N, label="t = $t", markersize=10, color=:red);

display(f2)