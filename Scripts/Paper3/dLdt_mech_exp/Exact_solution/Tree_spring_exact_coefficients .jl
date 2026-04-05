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

# Case: l₁(0) = l₂(0) = l₀
f = Figure(size=(1200,600))
ax = Axis(f[1,1], aspect=1, xlabel=L"\ell(0)", ylabel=L"\text{Coefficient values}", title=L"\text{Exact solution for } x_{1}(t)", limits=(-0.5,10.5, -5, 5))
# shading 
points = [Point2f(5, -11), Point2f(5, 11), Point2f(-2, 11), Point2f(-2, -11)]
poly!(ax, points, color = :green, alpha = 0.1)
points = [Point2f(5, -11), Point2f(5, 11), Point2f(12, 11), Point2f(12, -11)]
poly!(ax, points, color = :red, alpha = 0.1)
vlines!(ax, [5], linewidth=3, color=:black, linestyle=:dash, label=L"a")
lines!(ax, l₀, C₁(a, λ₁, λ₂, l₀, l₀), linewidth=3, label=L"C_{1}")
lines!(ax, l₀, C₂(a, λ₁, λ₂, l₀, l₀), linewidth=3, label=L"C_{2}")
axislegend(ax,position=:lt)

ax2 = Axis(f[1,2], aspect=1, xlabel=L"\ell(0)", title=L"\text{Exact solution for } x_{2}(t)", limits=(-0.5,10.5, -5, 5))
# shading 
points = [Point2f(5, -11), Point2f(5, 11), Point2f(-2, 11), Point2f(-2, -11)]
poly!(ax2, points, color = :green, alpha = 0.1)
points = [Point2f(5, -11), Point2f(5, 11), Point2f(12, 11), Point2f(12, -11)]
poly!(ax2, points, color = :red, alpha = 0.1)
vlines!(ax2, [5], linewidth=3, color=:black, linestyle=:dash)
lines!(ax2, l₀, (1 + λ₁) .* C₁(a, λ₁, λ₂, l₀, l₀), linewidth=3, label=L"(1+\lambda_{1}) C_{1}")
lines!(ax2, l₀, (1 + λ₂) .* C₂(a, λ₁, λ₂, l₀, l₀), linewidth=3, label=L"(1+\lambda_{2}) C_{2}")
axislegend(ax2,position=:lt)

display(f)
save("exact_sol_coefficients_uniform.png", f)


import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK


function compare_exact_discrete(k, a, η, l₁, l₂, l₃)

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
    if l₁ > a || l₂ > a
        axislegend(ax,position=:rt)
    else
        axislegend(ax,position=:rb)
    end
    return F
end

k = 2; a = 5.0; η = 1; l₁ = 2.0; l₂ = 2.0; l₃ = 5.0;
F = compare_exact_discrete(k, a, η, l₁, l₂, l₃)
save("Exact_vs_discrete_sol_$l₁"*"_$l₂.png", F)
