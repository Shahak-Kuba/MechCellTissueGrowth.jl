import MechCellTissueGrowth as MCTG

# Parameters
k = 5.0; η = 1.0; α = k/η; a = 5.0; L0 = 10.0; N = 101;

# Force and Diffusivity functions (Hookean spring model)
Ffunc = ρ -> k * (1/ρ - a)
Dfunc = ρ -> α / ρ^2 

# Force and Diffusivity functions (Nonlinear spring model)
#Ffunc = ρ -> k * (1/a - ρ)
#Dfunc = ρ -> α 

p = MCTG.FBParams(α=α, η=η, k=k, a=a, L0=L0, N=N, rBC=:free, lBC=:fixed, F=Ffunc, D=Dfunc)

q0_func = z -> 2.0 #* cos(z*π/2) + 1.0 #2.0
y0 = MCTG.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)

tspan = (0.0, 40.0)

# solve PDE
prob = ODEProblem(rhs!, y0, tspan, p)
sol = solve(prob, Rodas5P(), saveat=vcat([0.0:0.0001:40.0]...))