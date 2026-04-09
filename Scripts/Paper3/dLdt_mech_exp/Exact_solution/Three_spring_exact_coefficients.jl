# Eigenvalues: roots of λ³ + 5λ² + 6λ + 1 = 0  (via trigonometric method for 3 real roots)
θ  = acos(-1 / (2√7)) / 3
λ₁ = (2√7 / 3) * cos(θ)         - 5/3
λ₂ = (2√7 / 3) * cos(θ - 2π/3)  - 5/3
λ₃ = (2√7 / 3) * cos(θ - 4π/3)  - 5/3

C₁ = (a, λ₁, λ₂, λ₃, l₁, l₂, l₃) ->
    (2 .+ λ₁) .* ((l₁ .- a).*(1 .+ λ₂).*(1 .+ λ₃) .+ (l₂ .- a) .- (l₃ .- a).*(2 .+ λ₂).*(2 .+ λ₃)) ./
    ((λ₂ - λ₁) .* (λ₃ - λ₁))

C₂ = (a, λ₁, λ₂, λ₃, l₁, l₂, l₃) ->
    (2 .+ λ₂) .* ((l₁ .- a).*(1 .+ λ₁).*(1 .+ λ₃) .+ (l₂ .- a) .- (l₃ .- a).*(2 .+ λ₁).*(2 .+ λ₃)) ./
    ((λ₁ - λ₂) .* (λ₃ - λ₂))

C₃ = (a, λ₁, λ₂, λ₃, l₁, l₂, l₃) ->
    (2 .+ λ₃) .* ((l₁ .- a).*(1 .+ λ₁).*(1 .+ λ₂) .+ (l₂ .- a) .- (l₃ .- a).*(2 .+ λ₁).*(2 .+ λ₂)) ./
    ((λ₁ - λ₃) .* (λ₂ - λ₃))



# Solution components  (eigenvector i: [1, 1+λᵢ, (1+λᵢ)/(2+λᵢ)])
L₁ = (k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃) ->
    a .+ C₁(a,λ₁,λ₂,λ₃,l₁,l₂,l₃) .* exp.((k.*λ₁.*t)./η) .+
         C₂(a,λ₁,λ₂,λ₃,l₁,l₂,l₃) .* exp.((k.*λ₂.*t)./η) .+
         C₃(a,λ₁,λ₂,λ₃,l₁,l₂,l₃) .* exp.((k.*λ₃.*t)./η)

L₂ = (k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃) ->
    a .+ (1 .+ λ₁).*C₁(a,λ₁,λ₂,λ₃,l₁,l₂,l₃).*exp.((k.*λ₁.*t)./η) .+
         (1 .+ λ₂).*C₂(a,λ₁,λ₂,λ₃,l₁,l₂,l₃).*exp.((k.*λ₂.*t)./η) .+
         (1 .+ λ₃).*C₃(a,λ₁,λ₂,λ₃,l₁,l₂,l₃).*exp.((k.*λ₃.*t)./η)

L₃ = (k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃) ->
    a .+ ((1 .+ λ₁)./(2 .+ λ₁)).*C₁(a,λ₁,λ₂,λ₃,l₁,l₂,l₃).*exp.((k.*λ₁.*t)./η) .+
         ((1 .+ λ₂)./(2 .+ λ₂)).*C₂(a,λ₁,λ₂,λ₃,l₁,l₂,l₃).*exp.((k.*λ₂.*t)./η) .+
         ((1 .+ λ₃)./(2 .+ λ₃)).*C₃(a,λ₁,λ₂,λ₃,l₁,l₂,l₃).*exp.((k.*λ₃.*t)./η)

# uniform distribution l₁(0) = l₂(0)
l₀ = LinRange(0, 10, 101)

using CairoMakie
### Coefficient plots
a = 5.0;
# Case: l₁(0) = l₂(0) = l₀
f = Figure(size=(800,600))
ax = Axis(f[1,1], aspect=1, xlabel=L"\ell(0)", title=L"\text{Exact solution for } \ell_{3}(t)", limits=(-0.5,10.5, -5, 5))
# shading 
points = [Point2f(5, -11), Point2f(5, 11), Point2f(-2, 11), Point2f(-2, -11)]
poly!(ax, points, color = :green, alpha = 0.1)
points = [Point2f(5, -11), Point2f(5, 11), Point2f(12, 11), Point2f(12, -11)]
poly!(ax, points, color = :red, alpha = 0.1)
vlines!(ax, [5], linewidth=3, color=:black, linestyle=:dash)
lines!(ax, l₀, (1 + λ₁)/(2 + λ₁) .* C₁(a, λ₁, λ₂, λ₃, l₀, l₀, l₀), linewidth=3, label=L"\frac{1+\lambda_{1}}{2+\lambda_{1}} C_{1}")
lines!(ax, l₀, (1 + λ₂)/(2 + λ₂) .* C₂(a, λ₁, λ₂, λ₃, l₀, l₀, l₀), linewidth=3, label=L"\frac{1+\lambda_{2}}{2+\lambda_{2}} C_{2}")
lines!(ax, l₀, (1 + λ₃)/(2 + λ₃) .* C₃(a, λ₁, λ₂, λ₃, l₀, l₀, l₀), linewidth=3, label=L"\frac{1+\lambda_{3}}{2+\lambda_{3}} C_{3}")
Legend(f[1,2],ax)

display(f)
save("Three_springs_exact_sol_coefficients_uniform.png", f)


import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK


function compare_exact_discrete(k, a, η, l₁, l₂, l₃, legend_on=true)

    # solving discrete model
    Domain = MCTG.DomainProperties_t(N=3, m = 1, domain_type = "1D")
    CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t()
    Embed = MCTG.CellEvent_t()
    ProlifEmbed = MCTG.CellEvent_t()
    SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)
    IC = [l₁, l₁ + l₂, l₁ + l₂ + l₃]
    NumSaveTimePoints = 1001
    sol = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
    # discrete model cell length
    CL_1 = [1 ./ (sum(sol.Density[ii][1:Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_2 = [1 ./ (sum(sol.Density[ii][Domain.m + 1:2*Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_3 = [1 ./ (sum(sol.Density[ii][2*Domain.m + 1:end])/Domain.m) for ii in eachindex(sol.Density)]

    t = LinRange(0, SimTime.Tmax, 1001);
    y_min = minimum([l₁, l₂, l₃, a]) - 0.5
    y_max = maximum([l₁, l₂, l₃, a]) + 0.5
    F = Figure(size=(600,600))
    ax = Axis(F[1,1], aspect=1, xlabel=L"t", ylabel=L"\ell(t)", limits=(-0.5,10.5, y_min, y_max))
    lines!(ax, t, L₁(k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃), linewidth=3, label=L"\ell_{1}(t)", color=:navy)
    lines!(ax, t, L₂(k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃), linewidth=3, label=L"\ell_{2}(t)", color=:purple)
    lines!(ax, t, L₃(k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃), linewidth=3, label=L"\ell_{3}(t)", color=:darkorange)
    lines!(ax, sol.t, CL_1, linewidth = 3, linestyle = :dash, color = :cyan, label=L"\text{Discrete }\ell_{1}")
    lines!(ax, sol.t, CL_2, linewidth = 3, linestyle = :dash, color = :magenta, label=L"\text{Discrete }\ell_{2}")
    lines!(ax, sol.t, CL_3, linewidth = 3, linestyle = :dash, color = :red, label=L"\text{Discrete }\ell_{3}")
    if legend_on
        if l₁ > a || l₂ > a || l₃ > a
            axislegend(ax,position=:rt)
        else
            axislegend(ax,position=:rb)
        end
    end
    return F
end

function compare_exact_discrete(k, a, η, l₁, l₂, l₃, legend_on=true, title_txt="test")

    # solving discrete model
    Domain = MCTG.DomainProperties_t(N=3, m = 1, domain_type = "1D")
    CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t()
    Embed = MCTG.CellEvent_t()
    ProlifEmbed = MCTG.CellEvent_t()
    SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)
    IC = [l₁, l₁ + l₂, l₁ + l₂ + l₃]
    NumSaveTimePoints = 1001
    sol = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
    # discrete model cell length
    CL_1 = [1 ./ (sum(sol.Density[ii][1:Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_2 = [1 ./ (sum(sol.Density[ii][Domain.m + 1:2*Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_3 = [1 ./ (sum(sol.Density[ii][2*Domain.m + 1:end])/Domain.m) for ii in eachindex(sol.Density)]

    t = LinRange(0, SimTime.Tmax, 1001);
    y_min = minimum([l₁, l₂, l₃, a]) - 0.5
    y_max = maximum([l₁, l₂, l₃, a]) + 0.5
    F = Figure(size=(600,600))
    ax = Axis(F[1,1], aspect=1, xlabel=L"t", ylabel=L"\ell(t)", title=title_txt, limits=(-0.5,10.5, y_min, y_max))
    lines!(ax, t, L₁(k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃), linewidth=3, label=L"\ell_{1}(t)", color=:navy)
    lines!(ax, t, L₂(k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃), linewidth=3, label=L"\ell_{2}(t)", color=:purple)
    lines!(ax, t, L₃(k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃), linewidth=3, label=L"\ell_{3}(t)", color=:darkorange)
    lines!(ax, sol.t, CL_1, linewidth = 3, linestyle = :dash, color = :cyan, label=L"\text{Discrete }\ell_{1}")
    lines!(ax, sol.t, CL_2, linewidth = 3, linestyle = :dash, color = :magenta, label=L"\text{Discrete }\ell_{2}")
    lines!(ax, sol.t, CL_3, linewidth = 3, linestyle = :dash, color = :red, label=L"\text{Discrete }\ell_{3}")
    if legend_on
        if l₁ > a || l₂ > a || l₃ > a
            axislegend(ax,position=:rt)
        else
            axislegend(ax,position=:rb)
        end
    end
    return F
end

function compare_exact_discrete_include_IC(k, a, η, l₁, l₂, l₃, legend_on=true, title_txt="test")

    # solving discrete model
    Domain = MCTG.DomainProperties_t(N=3, m = 1, domain_type = "1D")
    CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t()
    Embed = MCTG.CellEvent_t()
    ProlifEmbed = MCTG.CellEvent_t()
    SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)
    IC = [l₁, l₁ + l₂, l₁ + l₂ + l₃]
    NumSaveTimePoints = 1001
    sol = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
    # discrete model cell length
    CL_1 = [1 ./ (sum(sol.Density[ii][1:Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_2 = [1 ./ (sum(sol.Density[ii][Domain.m + 1:2*Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_3 = [1 ./ (sum(sol.Density[ii][2*Domain.m + 1:end])/Domain.m) for ii in eachindex(sol.Density)]

    t = LinRange(0, SimTime.Tmax, 1001);
    y_min = minimum([l₁, l₂, l₃, a]) - 0.5
    y_max = maximum([l₁, l₂, l₃, a]) + 0.5
    F = Figure(size=(600,1200))
    ax = Axis(F[1,1], xlabel=L"x(t)", title="Initial condition", limits=(-0.5,IC[end] + 0.5, -0.1, 0.1))
    lines!(ax, [0; IC], [0, 0, 0, 0], color=:black, linewidth=3, label="IC")
    scatter!(ax, [0; IC], [0, 0, 0, 0], color=:red, markersize=20, label="IC")

    ax = Axis(F[2,1], xlabel=L"t", ylabel=L"\ell(t)", title=title_txt, limits=(-0.5,10.5, y_min, y_max))
    lines!(ax, t, L₁(k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃), linewidth=3, label=L"\ell_{1}(t)", color=:navy)
    lines!(ax, t, L₂(k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃), linewidth=3, label=L"\ell_{2}(t)", color=:purple)
    lines!(ax, t, L₃(k, η, t, a, λ₁, λ₂, λ₃, l₁, l₂, l₃), linewidth=3, label=L"\ell_{3}(t)", color=:darkorange)
    lines!(ax, sol.t, CL_1, linewidth = 3, linestyle = :dash, color = :cyan, label=L"\text{Discrete }\ell_{1}")
    lines!(ax, sol.t, CL_2, linewidth = 3, linestyle = :dash, color = :magenta, label=L"\text{Discrete }\ell_{2}")
    lines!(ax, sol.t, CL_3, linewidth = 3, linestyle = :dash, color = :red, label=L"\text{Discrete }\ell_{3}")
    if legend_on
        if l₁ > a || l₂ > a || l₃ > a
            axislegend(ax,position=:rt)
        else
            axislegend(ax,position=:rb)
        end
    end
    return F
end

k = 2; a = 5.0; η = 1; l₁ = 5.0; l₂ = 5.0; l₃ = 3.0;
F = compare_exact_discrete(k, a, η, l₁, l₂, l₃, true)
save("Three_springs_Exact_vs_discrete_sol_$l₁"*"_$l₂"*"_$l₃.png", F)

# special cases where c₁ ≠ 0, c₂ = c₃ = 0:
l₁ = 6.0; α = l₁ - a;
l₂ = a + (1 + λ₁) * α;
l₃ = a + (1 + λ₁) / (2 + λ₁) * α;
println("C₁ = ", C₁(a, λ₁, λ₂, λ₃, l₁, l₂, l₃), " C₂ = ", C₂(a, λ₁, λ₂, λ₃, l₁, l₂, l₃), " C₃ = ", C₃(a, λ₁, λ₂, λ₃, l₁, l₂, l₃))
F = compare_exact_discrete_include_IC(k, a, η, l₁, l₂, l₃, true, "C₁ ≠ 0")
save("Three_springs_E_vs_D_sol_c1_nonzero_tension.png", F)

# special cases where c₂ ≠ 0, c₁ = c₃ = 0:
l₁ = 6.0; α = l₁ - a;
l₂ = a + (1 + λ₂) * α;
l₃ = a + (1 + λ₂) / (2 + λ₂) * α;
println("C₁ = ", C₁(a, λ₁, λ₂, λ₃, l₁, l₂, l₃), " C₂ = ", C₂(a, λ₁, λ₂, λ₃, l₁, l₂, l₃), " C₃ = ", C₃(a, λ₁, λ₂, λ₃, l₁, l₂, l₃))
F = compare_exact_discrete_include_IC(k, a, η, l₁, l₂, l₃, false, "C₂ ≠ 0")
save("Three_springs_E_vs_D_sol_c2_nonzero_tension.png", F)

# special cases where c₃ ≠ 0, c₁ = c₂ = 0:
l₁ = 6.0; α = l₁ - a;
l₂ = a + (1 + λ₃) * α;
l₃ = a + (1 + λ₃) / (2 + λ₃) * α;
println("C₁ = ", C₁(a, λ₁, λ₂, λ₃, l₁, l₂, l₃), " C₂ = ", C₂(a, λ₁, λ₂, λ₃, l₁, l₂, l₃), " C₃ = ", C₃(a, λ₁, λ₂, λ₃, l₁, l₂, l₃))
F = compare_exact_discrete_include_IC(k, a, η, l₁, l₂, l₃, false, "C₃ ≠ 0")
save("Three_springs_E_vs_D_sol_c3_nonzero_tension.png", F)

# continuum simulation
# Parameters
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

# plotting initial condition
F = Figure(size=(600,600))
ax = Axis(F[1,1], xlabel=L"x", ylabel=L"\rho(x,0)", title="Initial condition")
lines!(ax, x, y0[1:end-1], color=:black, linewidth=3, label="IC")
scatter!(ax, x, y0[1:end-1], color=:red, markersize=20, label="IC")
display(F)
tspan = (0.0, 10.0)

# solve PDE
prob = ODEProblem(MCTG.rhs!, y0, tspan, p)
sol = solve(prob, Rodas5P(), saveat=vcat([0.0:0.001:10.0]...))

# plotting PDE solution
F = Figure(size=(600,600))
ax = Axis(F[1,1], xlabel=L"t", ylabel=L"L(t)", title="PDE solution")
lines!(ax, sol.t, [sol.u[ii][end] for ii in eachindex(sol.u)], color=:blue, linewidth=3, label="PDE")
display(F)

T = 2
F = Figure(size=(600,600))
ax = Axis(F[1,1], xlabel=L"x", ylabel=L"\rho(x,T)", title="T = $(sol.t[T])")
lines!(ax, collect(LinRange(0, 1.0, p.N)) .* sol.u[T][end], sol.u[T][1:end-1], color=:blue, linewidth=3, label="PDE")
display(F)

# Create animation of density evolution
fig = Figure(size=(600, 600))
ax_springs = Axis(fig[1, 1], xlabel=L"x", limits=(0, 16, -0.1, 0.1))
ax_density = Axis(fig[1, 1], xlabel=L"x", ylabel=L"\rho(x,t)", limits=(0, 16, 0, 0.4))

record(fig, "density_evolution.gif", 1:10:5001; framerate=30) do frame
    t_idx = frame
    t_val = sol.t[t_idx]
    # Discrete spring model
    empty!(ax_springs)

    # Density evolution
    empty!(ax)
    ax.title = "t = $(round(t_val, digits=3))"
    #ax.limits = (0, sol.u[end][end], 0, maximum(sol.u[1][1:end-1]) * 1.1)
    lines!(ax_density, collect(LinRange(0, 1.0, p.N)) .* sol.u[t_idx][end], sol.u[t_idx][1:end-1], 
           color=:blue, linewidth=3, label="PDE")
    axislegend(ax)

    # Evolution of spring lengths

end
