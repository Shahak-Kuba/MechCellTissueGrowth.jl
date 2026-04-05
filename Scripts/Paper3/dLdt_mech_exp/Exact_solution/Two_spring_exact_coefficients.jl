#C₁ = (a, λ₁, λ₂, l₁, l₂) -> ( (l₂ .- a) .- λ₂.*(l₁ .- a) ) ./ (√5)
#C₂ = (a, λ₁, λ₂, l₁, l₂) -> ( λ₁.*(l₁ .- a) .- (l₂ .- a) ) ./ (√5)

C₁ = (a, λ₁, λ₂, l₁, l₂) -> ( (l₁ .- a) .* (1 + λ₂) .- (l₂ .- a) ) ./ (λ₂ - λ₁)
C₂ = (a, λ₁, λ₂, l₁, l₂) -> ( (l₂ .- a) .- (l₁ .- a).*(1 + λ₁) ) ./ (λ₂ - λ₁)

a = 5;
λ₁ = (-3 + √5) / 2;
λ₂ = (-3 - √5) / 2;

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


a_diff = 2
# Case: l₁(0) < a
f2 = Figure(size=(1800,1300))
ax = Axis(f2[1,1], aspect=1, xlabel=L"\ell_{2}(0)", ylabel="Coefficients", title=L"\ell_{1}(0) = a - \delta a", limits=(-1,11, -10, 10), xticklabelsvisible=false)
lines!(ax, l₀, (1 + λ₁).*C₁(a, λ₁, λ₂, (a-a_diff).*ones(size(l₀)), l₀), linewidth=3, label=L"C_{1}")
lines!(ax, l₀, (1 + λ₂).*C₂(a, λ₁, λ₂, (a-a_diff).*ones(size(l₀)), l₀), linewidth=3, label=L"C_{2}")
axislegend(ax,position=:lt)
# Case: l₁(0) = a
ax = Axis(f2[1,2], aspect=1, xlabel=L"\ell_{2}(0)", title=L"\ell_{1}(0) = a", limits=(-1,11, -10, 10), xticklabelsvisible=false, yticklabelsvisible=false)
lines!(ax, l₀, (1 + λ₁).*C₁(a, λ₁, λ₂, a.*ones(size(l₀)), l₀), linewidth=3, label=L"C_{1}")
lines!(ax, l₀, (1 + λ₂).*C₂(a, λ₁, λ₂, a.*ones(size(l₀)), l₀), linewidth=3, label=L"C_{2}")
# Case: l₁(0) > a
ax = Axis(f2[1,3], aspect=1, xlabel=L"\ell_{2}(0)", title=L"\ell_{1}(0) = a + \delta a", limits=(-1,11, -10, 10), xticklabelsvisible=false, yticklabelsvisible=false)
lines!(ax, l₀, (1 + λ₁).*C₁(a, λ₁, λ₂, (a+a_diff).*ones(size(l₀)), l₀), linewidth=3, label=L"C_{1}")
lines!(ax, l₀, (1 + λ₂).*C₂(a, λ₁, λ₂, (a+a_diff).*ones(size(l₀)), l₀), linewidth=3, label=L"C_{2}")
# Case: l₂(0) < a
ax = Axis(f2[2,1], aspect=1, xlabel=L"\ell_{1}(0)", ylabel="Coefficients", title=L"\ell_{2}(0) = a - \delta a", limits=(-1,11, -10, 10))
lines!(ax, l₀, (1 + λ₁).*C₁(a, λ₁, λ₂, l₀, (a-a_diff).*ones(size(l₀))), linewidth=3, label=L"C_{1}")
lines!(ax, l₀, (1 + λ₂).*C₂(a, λ₁, λ₂, l₀, (a-a_diff).*ones(size(l₀))), linewidth=3, label=L"C_{2}")
# Case: l₂(0) = a
ax = Axis(f2[2,2], aspect=1, xlabel=L"\ell_{1}(0)", title=L"\ell_{2}(0) = a", limits=(-1,11, -10, 10), yticklabelsvisible=false)
lines!(ax, l₀, (1 + λ₁).*C₁(a, λ₁, λ₂, l₀, a.*ones(size(l₀))), linewidth=3, label=L"C_{1}")
lines!(ax, l₀, (1 + λ₂).*C₂(a, λ₁, λ₂, l₀, a.*ones(size(l₀))), linewidth=3, label=L"C_{2}")
# Case: l₂(0) > a
ax = Axis(f2[2,3], aspect=1, xlabel=L"\ell_{1}(0)", title=L"\ell_{2}(0) = a + \delta a", limits=(-1,11, -10, 10), yticklabelsvisible=false)
lines!(ax, l₀, (1 + λ₁).*C₁(a, λ₁, λ₂, l₀, (a+a_diff).*ones(size(l₀))), linewidth=3, label=L"C_{1}")
lines!(ax, l₀, (1 + λ₂).*C₂(a, λ₁, λ₂, l₀, (a+a_diff).*ones(size(l₀))), linewidth=3, label=L"C_{2}")

f3 = Figure(size=(1200,600))
δa_array = [-4, -2, -1, 1, 2, 4]
ax = Axis(f3[1,1], aspect=1, xlabel=L"\ell_{2}(0)", ylabel=L"C_{1}", title=L"\ell_{1}(0) - \delta a", limits=(-1,11, -10, 10))
for δa in δa_array
    lines!(ax, l₀, (1 + λ₁).*C₁(a, λ₁, λ₂, (a-δa).*ones(size(l₀)), l₀), linewidth=3, label="δa = $δa")
end
ax = Axis(f3[1,2], aspect=1, xlabel=L"\ell_{2}(0)", ylabel=L"C_{2}", title=L"\ell_{1}(0) - \delta a", limits=(-1,11, -10, 10))
for δa in δa_array
    lines!(ax, l₀, (1 + λ₂).*C₂(a, λ₁, λ₂, (a-δa).*ones(size(l₀)), l₀), linewidth=3, label="δa = $δa")
end

f4 = Figure(size=(1200,600))
ax = Axis(f4[1,1], aspect=1, xlabel=L"\ell_{1}(0)", ylabel=L"C_{1}", title=L"\ell_{2}(0) - \delta a", limits=(-1,11, -10, 10))
for δa in δa_array
    lines!(ax, l₀, (1 + λ₁).*C₁(a, λ₁, λ₂, l₀, (a-δa).*ones(size(l₀))), linewidth=3, label="δa = $δa")
end
ax = Axis(f4[1,2], aspect=1, xlabel=L"\ell_{1}(0)", ylabel=L"C_{2}", title=L"\ell_{2}(0) - \delta a", limits=(-1,11, -10, 10))
for δa in δa_array
    lines!(ax, l₀, (1 + λ₂).*C₂(a, λ₁, λ₂, l₀, (a-δa).*ones(size(l₀))), linewidth=3, label="δa = $δa")
end

display(f)
display(f2)
display(f3)
display(f4)

save("exact_sol_coefficients_uniform.png", f)
save("exact_sol_coefficients_vary.png", f2)


import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK


function compare_exact_discrete(k, a, η, l₁, l₂)
    # Exact solutions
    L₁ = (k, η, t, a, λ₁, λ₂, l₁, l₂) -> a .+ C₁(a, λ₁, λ₂, l₁, l₂) .* exp.((k.*λ₁.*t)./η) .+ C₂(a, λ₁, λ₂, l₁, l₂) .* exp.((k.*λ₂.*t)./η);
    L₂ = (k, η, t, a, λ₁, λ₂, l₁, l₂) -> a .+ (1 .+ λ₁).*C₁(a, λ₁, λ₂, l₁, l₂).*exp.((k.*λ₁.*t)./η) .+ (1 .+ λ₂).*C₂(a, λ₁, λ₂, l₁, l₂) .* exp.((k.*λ₂.*t)./η);

    # solving discrete model
    Domain = MCTG.DomainProperties_t(N=2, m = 1, domain_type = "1D")
    CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t()
    Embed = MCTG.CellEvent_t()
    ProlifEmbed = MCTG.CellEvent_t()
    SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)
    IC = [l₁, l₁ + l₂]
    NumSaveTimePoints = 1001
    sol = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
    # discrete model cell length
    CL_1 = [1 ./ (sum(sol.Density[ii][1:Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_2 = [1 ./ (sum(sol.Density[ii][Domain.m + 1:end])/Domain.m) for ii in eachindex(sol.Density)]

    t = LinRange(0, SimTime.Tmax, 1001);
    y_min = minimum([l₁, l₂, a]) - 0.5
    y_max = maximum([l₁, l₂, a]) + 0.5
    F = Figure(size=(600,600))
    ax = Axis(F[1,1], aspect=1, xlabel=L"t", ylabel=L"\ell(t)", limits=(-0.5,10.5, y_min, y_max))
    lines!(ax, t, L₁(k, η, t, a, λ₁, λ₂, l₁, l₂), linewidth=3, label=L"\ell_{1}(t)", color=:navy)
    lines!(ax, t, L₂(k, η, t, a, λ₁, λ₂, l₁, l₂), linewidth=3, label=L"\ell_{2}(t)", color=:darkorange)
    lines!(ax, sol.t, CL_1, linewidth = 3, linestyle = :dash, color = :cyan, label=L"\text{Discrete }\ell_{1}")
    lines!(ax, sol.t, CL_2, linewidth = 3, linestyle = :dash, color = :red, label=L"\text{Discrete }\ell_{2}")
    if l₁ > a || l₂ > a
        axislegend(ax,position=:rt)
    else
        axislegend(ax,position=:rb)
    end
    return F
end

k = 2; a = 5.0; η = 1; l₁ = 8.0; l₂ = 5.0;
F = compare_exact_discrete(k, a, η, l₁, l₂)
save("Exact_vs_discrete_sol_$l₁"*"_$l₂.png", F)
