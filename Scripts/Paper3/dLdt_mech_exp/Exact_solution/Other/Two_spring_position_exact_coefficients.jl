C₁ = (a, λ₁, λ₂, x₁, x₂) -> ( (x₂ .- (2 .* a)) .- (2 .+ λ₂) .*(x₁ .- a)) ./ (λ₁ - λ₂)
C₂ = (a, λ₁, λ₂, x₁, x₂) -> ( (x₁ .- a).*(2 .+ λ₁) .- ( x₂ .- 2 .* a) ) ./ (λ₁ .- λ₂)

a = 5;
λ₁ = (-3 + √5) / 2;
λ₂ = (-3 - √5) / 2;

# uniform distribution l₁(0) = l₂(0)
l0 = LinRange(1, 9, 101)
x₁0 = l0;
x₂0 = 2 .* l0;


using CairoMakie
### Coefficient plots

# Case: l₁(0) = l₂(0) = l₀
f = Figure(size=(600,600))
ax = Axis(f[1,1], aspect=1, xlabel=L"\ell(0)", ylabel="Coefficients", title=L"\ell(0) = \ell_{1}(0) = \ell_{2}(0)", limits=(-1,11, -10, 10))
# shading 
points = [Point2f(5, -11), Point2f(5, 11), Point2f(-2, 11), Point2f(-2, -11)]
poly!(ax, points, color = :black, alpha = 0.05)
vlines!(ax, [5], linewidth=3, color=:black, linestyle=:dash, label=L"a")
lines!(ax, x₁0, C₁(a, λ₁, λ₂, x₁0, x₁0), linewidth=3, label=L"C_{1}")
lines!(ax, x₁0, C₂(a, λ₁, λ₂, x₁0, x₁0), linewidth=3, label=L"C_{2}")
axislegend(ax,position=:lt)

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
#display(f3)

# Exact solutions
L₁ = (k, η, t, a, λ₁, λ₂, x₁, x₂) -> a .+ C₁(a, λ₁, λ₂, x₁, x₂) .* exp.((k.*λ₁.*t)./η) .+ C₂(a, λ₁, λ₂, x₁, x₂) .* exp.((k.*λ₂.*t)./η);
L₂ = (k, η, t, a, λ₁, λ₂, x₁, x₂) -> (2 .* a) .+ (2 .+ λ₁).*C₁(a, λ₁, λ₂, x₁, x₂).*exp.((k.*λ₁.*t)./η) .+ (2 .+ λ₂).*C₂(a, λ₁, λ₂, x₁, x₂) .* exp.((k.*λ₂.*t)./η);
k = 1.0; η = 1.0; x₁ = 4.0; x₂ = 8.0;


# Discrete model
import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK
# solving discrete model
Domain = MCTG.DomainProperties_t(N=2, m = 1, domain_type = "1D")
CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
Prolif = MCTG.CellEvent_t()
Death = MCTG.CellEvent_t()
Embed = MCTG.CellEvent_t()
ProlifEmbed = MCTG.CellEvent_t()
SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)

#IC = [x₁, x₂]
l0 = x₂
IC = reverse(collect(LinRange(0, l0, Domain.N*Domain.m + 1)))
pop!(IC)
reverse!(IC)
NumSaveTimePoints = 1001
sol = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
# discrete model cell length
CB_1 = [sol.u[ii][1,Domain.m + 1] for ii in eachindex(sol.t)]
CB_2 = [sol.u[ii][1,Domain.N * Domain.m + 1] for ii in eachindex(sol.t)]


t = LinRange(0, SimTime.Tmax, 1001);

F = Figure(size=(600,600))
ax = Axis(F[1,1], aspect=1, xlabel=L"t", ylabel=L"x\;(t)", limits=(-1,11, 1, 11))
lines!(ax, t, L₁(k, η, t, a, λ₁, λ₂, x₁, x₂), linewidth=3, label=L"\x_{1}(t)", color=:navy)
lines!(ax, t, L₂(k, η, t, a, λ₁, λ₂, x₁, x₂), linewidth=3, label=L"\x_{2}(t)", color=:darkorange)
#lines!(ax, t, L₁(k, η, t, a, λ₁, λ₂, l₁, l₂) .+ L₂(k, η, t, a, λ₁, λ₂, l₁, l₂), linewidth=3, label=L"\ell_{1}(t) + \ell_{2}(t)", color=:green)
lines!(ax, sol.t, CB_1, linewidth = 3, linestyle = :dash, color = :cyan)
lines!(ax, sol.t, CB_2, linewidth = 3, linestyle = :dash, color = :red)
#axislegend(ax,position=:rb)
display(F)

