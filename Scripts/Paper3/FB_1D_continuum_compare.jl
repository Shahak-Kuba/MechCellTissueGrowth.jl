import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq

# Parameters
k = 5.0; η = 1.0; α = k/η; a = 2.0; L0 = 10.0; N = 1001; Δx = 1/(N-1); L0 = 10.0; q₀ = 1.0; 
N_IC = Int(L0 * q₀);
rBC=:free; 
tspan = (0.0, 10.0);

# Force and Diffusivity functions (Hookean spring model)
Ffunc = ρ -> k * (1/ρ - a)
Dfunc = ρ -> α / ρ^2 
# Proliferation and apoptosis functions
Pfunc = ρ -> 0.0
Afunc = ρ -> 0.0
#initial condition function
q0_func = z -> q₀ #* cos(z*π/2) + 1.0 #2.0

# solve original Baker PDE without correction
p = MCTG.FBParams(α=α, η=η, k=k, a=a, L0=L0, N=N, rBC=:free, lBC=:fixed, F=Ffunc, D=Dfunc, P=Pfunc, A=Afunc)
y0 = MCTG.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)
# solve PDE
prob = ODEProblem(MCTG.rhs_Baker!, y0, tspan, p)
sol_Baker = solve(prob, Rodas5P(), saveat=vcat([0.0:0.0001:tspan[2]]...))

# solve PDE with improper correction
p = MCTG.FBParams(α=α, η=η, k=k, a=a, L0=L0, N=N, rBC=:free, lBC=:fixed, F=Ffunc, D=Dfunc, P=Pfunc, A=Afunc)
y0 = MCTG.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)
# solve PDE
prob = ODEProblem(MCTG.rhs!, y0, tspan, p)
sol_new = solve(prob, Rodas5P(), saveat=vcat([0.0:0.0001:tspan[2]]...))

# solve PDE with correction
p = MCTG.FBParams_w_correction(α, k, a, η, Δx, L0, N, N_IC, q₀, rBC, Ffunc, Dfunc, Pfunc, Afunc)
y0 = MCTG.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)
# solve PDE
prob = ODEProblem(MCTG.rhs_with_correction!, y0, tspan, p)
sol_new_correction = solve(prob, Rodas5P(), saveat=vcat([0.0:0.0001:tspan[2]]...))

# plotting restults
f = Figure(size=(650,500))
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"L(t)", title="Comparison of PDE solutions")
lines!(ax, sol_Baker.t, [sol_Baker.u[j][end] for j in 1:length(sol_Baker.t)], linewidth=3, color=:blue, label="Baker")
lines!(ax, sol_new.t, [sol_new.u[j][end] for j in 1:length(sol_new.t)], linewidth=3, color=:red, label="PDE with O(m) term")
lines!(ax, sol_new_correction.t, [sol_new_correction.u[j][end] for j in 1:length(sol_new_correction.t)], linewidth=3, linestyle=:dash, color=:black, label="PDE with O(1) sum")
axislegend(ax, position=:rb)
display(f)
save("PDE_comparison.png", f)

mF_q_term_sol_new = [(100/η) .* Ffunc(sol_new.u[j][end-1]) for j in 1:length(sol_new.t)]
mF_q_term_correction = [(2 / (N_IC * η)) .* Ffunc(q₀) .* sum(exp.(-(k .* ((2 .* s) - 1).^2 .* π.^2 .* sol_new_correction.t[j]) ./ (4 .* (N_IC.^2) .* η)) for s in 1:500) for j in 1:length(sol_new_correction.t)]
f = Figure(size=(650,500))
ax = Axis(f[1,1], xlabel=L"t", ylabel=L"-m/η * F(q)", title="Comparison of mF_q_term")
lines!(ax, sol_new.t[2:end], -(mF_q_term_sol_new-mF_q_term_correction)[2:end], color=:blue)
display(f)