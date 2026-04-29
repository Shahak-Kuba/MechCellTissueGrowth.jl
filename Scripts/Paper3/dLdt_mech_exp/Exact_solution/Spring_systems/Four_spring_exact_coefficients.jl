# Eigenvalues: λ₁ = -1 (exact), λ₂₋₄ from t³ - 3t - 1 = 0 via cos(3θ) = 1/2
λ₁ =  -1.0
λ₂ = 2cos(π/9)  - 2
λ₃ = 2cos(5π/9) - 2
λ₄ = 2cos(7π/9) - 2

# Eigenvector components (from forward recurrence on rows of A - λᵢI)
p₂ = λ ->  1 + λ
p₃ = λ ->  λ^2 + 3λ + 1
p₄ = λ ->  λ^3 + 5λ^2 + 6λ + 1

# Constants (using orthogonality of eigenvectors since A is symmetric: cᵢ = vᵢᵀα / ‖vᵢ‖²)
nsq = λ -> 1 + p₂(λ)^2 + p₃(λ)^2 + p₄(λ)^2

C = (a, λᵢ, l₁, l₂, l₃, l₄) ->
    ((l₁ .- a) .+ p₂(λᵢ).*(l₂ .- a) .+ p₃(λᵢ).*(l₃ .- a) .+ p₄(λᵢ).*(l₄ .- a)) ./ nsq(λᵢ)

# Solution components
L₁ = (k, η, t, a, λ₁, λ₂, λ₃, λ₄, l₁, l₂, l₃, l₄) ->
    a .+ C(a,λ₁,l₁,l₂,l₃,l₄).*exp.((k.*λ₁.*t)./η) .+
         C(a,λ₂,l₁,l₂,l₃,l₄).*exp.((k.*λ₂.*t)./η) .+
         C(a,λ₃,l₁,l₂,l₃,l₄).*exp.((k.*λ₃.*t)./η) .+
         C(a,λ₄,l₁,l₂,l₃,l₄).*exp.((k.*λ₄.*t)./η)

L₂ = (k, η, t, a, λ₁, λ₂, λ₃, λ₄, l₁, l₂, l₃, l₄) ->
    a .+ p₂(λ₁).*C(a,λ₁,l₁,l₂,l₃,l₄).*exp.((k.*λ₁.*t)./η) .+
         p₂(λ₂).*C(a,λ₂,l₁,l₂,l₃,l₄).*exp.((k.*λ₂.*t)./η) .+
         p₂(λ₃).*C(a,λ₃,l₁,l₂,l₃,l₄).*exp.((k.*λ₃.*t)./η) .+
         p₂(λ₄).*C(a,λ₄,l₁,l₂,l₃,l₄).*exp.((k.*λ₄.*t)./η)

L₃ = (k, η, t, a, λ₁, λ₂, λ₃, λ₄, l₁, l₂, l₃, l₄) ->
    a .+ p₃(λ₁).*C(a,λ₁,l₁,l₂,l₃,l₄).*exp.((k.*λ₁.*t)./η) .+
         p₃(λ₂).*C(a,λ₂,l₁,l₂,l₃,l₄).*exp.((k.*λ₂.*t)./η) .+
         p₃(λ₃).*C(a,λ₃,l₁,l₂,l₃,l₄).*exp.((k.*λ₃.*t)./η) .+
         p₃(λ₄).*C(a,λ₄,l₁,l₂,l₃,l₄).*exp.((k.*λ₄.*t)./η)

L₄ = (k, η, t, a, λ₁, λ₂, λ₃, λ₄, l₁, l₂, l₃, l₄) ->
    a .+ p₄(λ₁).*C(a,λ₁,l₁,l₂,l₃,l₄).*exp.((k.*λ₁.*t)./η) .+
         p₄(λ₂).*C(a,λ₂,l₁,l₂,l₃,l₄).*exp.((k.*λ₂.*t)./η) .+
         p₄(λ₃).*C(a,λ₃,l₁,l₂,l₃,l₄).*exp.((k.*λ₃.*t)./η) .+
         p₄(λ₄).*C(a,λ₄,l₁,l₂,l₃,l₄).*exp.((k.*λ₄.*t)./η)


# uniform distribution l₁(0) = l₂(0)
l₀ = LinRange(0, 10, 101)

using CairoMakie
import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK


function compare_exact_discrete(k, a, η, l₁, l₂, l₃, l₄, legend_on=true)

    # solving discrete model
    Domain = MCTG.DomainProperties_t(N=4, m = 1, domain_type = "1D")
    CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t()
    Embed = MCTG.CellEvent_t()
    ProlifEmbed = MCTG.CellEvent_t()
    SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)
    IC = [l₁, l₁ + l₂, l₁ + l₂ + l₃, l₁ + l₂ + l₃ + l₄]
    NumSaveTimePoints = 1001
    sol = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
    # discrete model cell length
    CL_1 = [1 ./ (sum(sol.Density[ii][1:Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_2 = [1 ./ (sum(sol.Density[ii][Domain.m + 1:2*Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_3 = [1 ./ (sum(sol.Density[ii][2*Domain.m + 1:3*Domain.m])/Domain.m) for ii in eachindex(sol.Density)]
    CL_4 = [1 ./ (sum(sol.Density[ii][3*Domain.m + 1:end])/Domain.m) for ii in eachindex(sol.Density)]


    t = LinRange(0, SimTime.Tmax, 1001);
    y_min = minimum([l₁, l₂, l₃, l₄, a]) - 0.5
    y_max = maximum([l₁, l₂, l₃, l₄, a]) + 0.5
    F = Figure(size=(600,600))
    ax = Axis(F[1,1], aspect=1, xlabel=L"t", ylabel=L"\ell(t)", limits=(-0.5,10.5, y_min, y_max))
    lines!(ax, t, L₁(k, η, t, a, λ₁, λ₂, λ₃, λ₄, l₁, l₂, l₃, l₄), linewidth=3, label=L"\ell_{1}(t)", color=:navy)
    lines!(ax, t, L₂(k, η, t, a, λ₁, λ₂, λ₃, λ₄, l₁, l₂, l₃, l₄), linewidth=3, label=L"\ell_{2}(t)", color=:purple)
    lines!(ax, t, L₃(k, η, t, a, λ₁, λ₂, λ₃, λ₄, l₁, l₂, l₃, l₄), linewidth=3, label=L"\ell_{3}(t)", color=:darkorange)
    lines!(ax, t, L₄(k, η, t, a, λ₁, λ₂, λ₃, λ₄, l₁, l₂, l₃, l₄), linewidth=3, label=L"\ell_{4}(t)", color=:green)
    lines!(ax, sol.t, CL_1, linewidth = 3, linestyle = :dash, color = :cyan, label=L"\text{Discrete }\ell_{1}")
    lines!(ax, sol.t, CL_2, linewidth = 3, linestyle = :dash, color = :magenta, label=L"\text{Discrete }\ell_{2}")
    lines!(ax, sol.t, CL_3, linewidth = 3, linestyle = :dash, color = :red, label=L"\text{Discrete }\ell_{3}")
    lines!(ax, sol.t, CL_4, linewidth = 3, linestyle = :dash, color = :black, label=L"\text{Discrete }\ell_{4}")
    if legend_on
        if l₁ > a || l₂ > a || l₃ > a || l₄ > a
            axislegend(ax,position=:rt)
        else
            axislegend(ax,position=:rb)
        end
    end
    return F
end

# Parameters
k = 2; a = 5.0; η = 1; l₁ = 6.0; l₂ = 6.0; l₃ = 6.0; l₄ = 6.0;
F = compare_exact_discrete(k, a, η, l₁, l₂, l₃, l₄, false)
save("Four_springs_Exact_vs_discrete_sol_$l₁"*"_$l₂"*"_$l₃"*"_$l₄.png", F)

