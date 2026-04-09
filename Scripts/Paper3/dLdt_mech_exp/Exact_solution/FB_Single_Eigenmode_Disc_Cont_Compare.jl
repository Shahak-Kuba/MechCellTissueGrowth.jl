import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK

Single_Eigenmode = "C1" # options: "C1", "C2", "C3"

l₂ = 0.0; l₃ = 0.0;

l₁ = 6.0; α = l₁ - a;
if Single_Eigenmode == "C1"
    # special cases where c₁ ≠ 0, c₂ = c₃ = 0:
    l₂ = a + (1 + λ₁) * α;
    l₃ = a + (1 + λ₁) / (2 + λ₁) * α;
elseif Single_Eigenmode == "C2"
    # special cases where c₂ ≠ 0, c₁ = c₃ = 0:
    l₂ = a + (1 + λ₂) * α;
    l₃ = a + (1 + λ₂) / (2 + λ₂) * α;
elseif Single_Eigenmode == "C3"
    # special cases where c₃ ≠ 0, c₁ = c₂ = 0:
    l₂ = a + (1 + λ₃) * α;
    l₃ = a + (1 + λ₃) / (2 + λ₃) * α;
else
    error("Single_Eigenmode variable not recognised")
end 

## solving discrete model
Domain = MCTG.DomainProperties_t(N=3, m = 1, domain_type = "1D")
CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
Prolif = MCTG.CellEvent_t()
Death = MCTG.CellEvent_t()
Embed = MCTG.CellEvent_t()
ProlifEmbed = MCTG.CellEvent_t()
SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)
IC = [l₁, l₁ + l₂, l₁ + l₂ + l₃]
NumSaveTimePoints = 1001
sol_disc = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);

## solving continuum model

α = k/η;  N = 1001;
# Force and Diffusivity functions (Hookean spring model)
Ffunc = ρ -> k * (1/ρ - a)
Dfunc = ρ -> α / ρ^2 

# Proliferation and apoptosis functions
Pfunc = ρ -> 0.0
Afunc = ρ -> 0.0

# making initial conditions from lengths
Cell_Lengths = [l₁, l₂, l₃]
L0 = cumsum(Cell_Lengths)[end]

p = MCTG.FBParams(α=α, η=η, k=k, a=a, L0=L0, N=N, rBC=:free, lBC=:fixed, F=Ffunc, D=Dfunc, P=Pfunc, A=Afunc)
# Initial conditions
Cell_densities = 1 ./ Cell_Lengths
Normalized_Cell_Lengths = cumsum(Cell_Lengths ./ L0)
x = LinRange(0, 1, p.N)
y0 = zeros(p.N + 1)
for (ii, norm_cell_length) in enumerate(Normalized_Cell_Lengths)
    for jj in eachindex(x)
        if y0[jj] < 0.1 && x[jj] <= norm_cell_length
            y0[jj] = Cell_densities[ii]
        end
    end
end
y0[end] = L0

# solve PDE
prob = ODEProblem(MCTG.rhs!, y0, tspan, p)
sol_cont = solve(prob, Rodas5P(), saveat=vcat([0.0:0.01:10.0]...))

## Plotting
# Create animation of density evolution
fig = Figure(size=(1200, 1300))
ax_springs = Axis(fig[1, 1:2], xlabel=L"x", limits=(0, L0 + 0.5, -0.1, 0.1))
ax_density = Axis(fig[2, 1], xlabel=L"x", ylabel=L"\rho(x,t)", limits=(0, L0 + 0.5, 0, 0.4))
ax_L_t = Axis(fig[2, 2], xlabel=L"t", ylabel=L"L(t)", limits=(0, 10, 14, L0 + 0.5))

# precacluting L(t) for discrete and continuum models
L_disc = [sol_disc.u[i][1, end] for i in eachindex(sol_disc.t)]
L_cont = [sol_cont.u[i][end] for i in eachindex(sol_cont.t)]

record(fig, "density_evolution_C1.gif", 1:501; framerate=30) do frame
    t_idx = frame
    t_val = sol.t[t_idx]
    # Discrete spring model
    empty!(ax_springs)
    x_springs = sol_disc.u[t_idx][1,:]
    y_springs = zeros(length(x_springs))
    lines!(ax_springs, x_springs, y_springs, linewidth=3, color=:black)
    scatter!(ax_springs, x_springs, y_springs, markersize=20, marker=:circle, color=:red)

    # Density evolution
    empty!(ax_density)
    ax.title = "t = $(round(t_val, digits=3))"
    disc_density = sol_disc.Density[t_idx]
    mid_disc_spring_pos = (x_springs[1:end-1] + x_springs[2:end]) ./ 2
    #ax.limits = (0, sol.u[end][end], 0, maximum(sol.u[1][1:end-1]) * 1.1)
    scatter!(ax_density, mid_disc_spring_pos, disc_density, markersize=20, marker=:circle, color=:red, label="Discrete")
    lines!(ax_density, collect(LinRange(0, 1.0, p.N)) .* sol.u[t_idx][end], sol.u[t_idx][1:end-1], 
           color=:blue, linewidth=3, label="Continuum")
    axislegend(ax_density)

    # Evolution of spring lengths
    empty!(ax_L_t)
    lines!(ax_L_t, sol_disc.t[1:t_idx], L_disc[1:t_idx], label="Discrete", linewidth=3, color=:red)
    lines!(ax_L_t, sol_cont.t[1:t_idx], L_cont[1:t_idx], label="Continuum", linewidth=3, color=:blue)
    axislegend(ax_L_t)
end