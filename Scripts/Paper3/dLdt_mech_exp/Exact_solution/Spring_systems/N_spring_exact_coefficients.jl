
θ  = (p, N) -> (2p - 1) * π / (2N + 1)
λ  = (p, N) -> 2cos(θ(p, N)) - 2
vN = (p, N) -> cos((N - 0.5) * θ(p, N)) / cos(0.5 * θ(p, N))

C  = (a, p, N, l₀) -> (4cos.(0.5 .* θ.(p, N)) ./ (2N .+ 1)) *
                        sum(cos.((j .- 0.5) .* θ.(p, N)) .* (l₀[j] .- a) for j in 1:N)

Lₙ = (k, η, t, a, N, l₀) ->
    a .+ sum(C(a, p, N, l₀) .* vN(p, N) .* exp.((k .* λ(p, N) .* t) ./ η) for p in 1:N)

λ_Kouachi = (p, N) -> -2 + 2cos((2*p*π) / (2*N+1))

using CairoMakie

# figure to compare eigenvalues 
f = Figure(size=(600,600))
ax = Axis(f[1,1], aspect=1, xlabel=L"p", ylabel=L"\lambda_{p}")
N_array = [10, 50, 100]
for N in N_array
    p_array = 1:N
    λ_array = [λ(p, N) for p in p_array]
    λ_K_array = [λ_Kouachi(p, N) for p in p_array]
    if N == N_array[1]
        legend_on = true
    else
        legend_on = false
    end
    scatter!(ax, p_array, λ_array, markersize=10,marker=:cross, color=:navy, label=L"\text{Derived}")
    scatter!(ax, p_array, λ_K_array, markersize=10, marker=:star5, color=:darkorange, label=L"\text{Kouachi.~2006}")
    if legend_on
        axislegend(ax, position=:rt)
    end
end
#axislegend(ax, position=:rt)
display(f)
save("eigenvalues_comparison.png", f)

import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using QuadGK


function compare_exact_discrete(N::Int64, k, a, η, l_array, legend_on=true)

    # solving discrete model
    Domain = MCTG.DomainProperties_t(N=N, m = 1, domain_type = "1D")
    CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t()
    Embed = MCTG.CellEvent_t()
    ProlifEmbed = MCTG.CellEvent_t()
    SimTime = MCTG.SimTime_t(Tmax=10, δt=0.001, event_δt=0.001)
    IC = cumsum(l_array)
    NumSaveTimePoints = 1001
    sol = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
    # discrete model cell length
    CL_N = [1 ./ (sum(sol.Density[ii][N-1*Domain.m + 1:end])/Domain.m) for ii in eachindex(sol.Density)]


    t = LinRange(0, SimTime.Tmax, 1001);
    y_min = minimum([l_array; a]) - 0.5
    y_max = maximum([l_array; a]) + 0.5
    F = Figure(size=(600,600))
    ax = Axis(F[1,1], aspect=1, xlabel=L"t", ylabel=L"\ell(t)", title ="N = $N", limits=(-0.5,10.5, y_min, y_max))
    lines!(ax, t, Lₙ(k, η, t, a, N, l_array), linewidth=3, label=L"\ell_{N}(t)", color=:navy)
    lines!(ax, sol.t, CL_N, linewidth = 3, linestyle = :dash, color = :cyan, label=L"\text{Discrete }\ell_{N}")
    if legend_on
        if sum(l_array .> a) < 0 
            axislegend(ax,position=:rt)
        else
            axislegend(ax,position=:rb)
        end
    end
    return F
end

function compare_exact_discrete(N::Vector{Int64}, k, a, η, l₀, legend_on=true)
    F = Figure(size=(600,600))
    y_min = minimum([l₀; a]) - 0.5
    y_max = maximum([l₀; a]) + 0.5
    ax = Axis(F[1,1], aspect=1, xlabel=L"t", ylabel=L"\ell_{N}(t)", limits=(-0.5,10.5, 2.0, y_max))
    for N_val in N
        # solving discrete model
        Domain = MCTG.DomainProperties_t(N=N_val, m = 10, domain_type = "1D")
        CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
        Prolif = MCTG.CellEvent_t()
        Death = MCTG.CellEvent_t()
        Embed = MCTG.CellEvent_t()
        ProlifEmbed = MCTG.CellEvent_t()
        SimTime = MCTG.SimTime_t(Tmax=10, δt=0.0001, event_δt=0.001)

        l_spring = l₀ / (Domain.N * Domain.m + 1) # length of each spring such that total length is l₀
        l_array = fill(l_spring, N_val)

        if Domain.m == 1
            IC = cumsum(l_array)
        else
            IC_CB = [0; cumsum(l_array)]
            IC = vcat([collect(LinRange(IC_CB[ii], IC_CB[ii+1], Domain.m+1))[2:end] for ii in 1:length(IC_CB)-1]...)
            #println("IC: ", IC)
        end
        NumSaveTimePoints = 1001
        sol = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);
        # discrete model cell length
        CL_N = [1 ./ (sum(sol.Density[ii][(N_val-1)*Domain.m + 1:end])/Domain.m) for ii in eachindex(sol.Density)]


        t = LinRange(0, SimTime.Tmax, 1001);
        exact_cell_length = zeros(size(Lₙ(k, η, t, a, N_val, l_array)))
        for ii in N_val-Domain.m+1:N_val
            exact_cell_length = exact_cell_length .+ Lₙ(k*Domain.m, η/Domain.m, t, a/Domain.m, ii, l_array)
        end
        #lines!(ax, t, Lₙ(k, η, t, a, N_val, l_array), linewidth=3, label="N = $N_val")
        lines!(ax, t, exact_cell_length, linewidth=3, label="N = $N_val")
        lines!(ax, sol.t, CL_N, linewidth = 3, linestyle = :dash, color=:black)
        if N_val == N[end]
            if legend_on
                if sum(l_array .> a) < 0 
                    axislegend(ax,position=:rt)
                else
                    axislegend(ax,position=:rb)
                end
            end
        end
    end
    return F
end
# Parameters
k = 3; a = 5.0; η = 1; N = 50;
N_array = [2, 5, 50]
l₀ = 8.0
F = compare_exact_discrete(N_array, k, a, η, l₀, true)
save("N_springs_Exact_vs_discrete_sol_$l₀.png", F)


# looking at eigenvalues vs coefficients for various initial conditions for N = 50
k = 3; a = 5.0; η = 1; N = 5000;
l₀_array = [2.0]

C_values = zeros(1,N)
λ_values = zeros(1,N)

f = Figure(size=(600,600))
ax = Axis(f[1,1], aspect=1, xlabel=L"\lambda_{p}", ylabel=L"C_{p}", title=L"\text{Coefficients for } \ell_{N}(t)")
for l₀ in l₀_array
    l_array = fill(l₀, N)
    C_values = [C(a, p, N, l_array) for p in 1:N]
    λ_values = [λ(p, N) for p in 1:N]
    scatter!(ax, λ_values, C_values, label = "l₀ = $l₀")
end
#axislegend(ax, position=:rt)
display(f)


count(>=(0.001), abs.(C_values))
"""
f = Figure(size=(600,600))
ax = Axis(f[1,1], aspect=1, xlabel=L"\lambda_{p}", ylabel=L"C_{p}", title=L"\text{Coefficients for } \ell_{N}(t)")
scatter!(ax, λ_values, C_values, markersize=10, color=:black)
display(f)
"""