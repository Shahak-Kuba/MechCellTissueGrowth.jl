import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq
using Printf
using QuadGK

Single_Eigenmode = "C3" # options: "C1", "C2", "C3"

function generate_3_spring_FB_animation(k, a, η, l₁, animation_name; maxT=10.0, special_type="None")
    
    θ  = acos(-1 / (2√7)) / 3
    λ₁ = (2√7 / 3) * cos(θ)         - 5/3
    λ₂ = (2√7 / 3) * cos(θ - 2π/3)  - 5/3
    λ₃ = (2√7 / 3) * cos(θ - 4π/3)  - 5/3

    α = l₁ - a;
    l₂ = 0.0; l₃ = 0.0;
    if special_type == "C1"
        # special cases where c₁ ≠ 0, c₂ = c₃ = 0:
        l₂ = a + (1 + λ₁) * α;
        l₃ = a + (1 + λ₁) / (2 + λ₁) * α;
    elseif special_type == "C2"
        # special cases where c₂ ≠ 0, c₁ = c₃ = 0:
        l₂ = a + (1 + λ₂) * α;
        l₃ = a + (1 + λ₂) / (2 + λ₂) * α;
    elseif special_type == "C3"
        # special cases where c₃ ≠ 0, c₁ = c₂ = 0:
        l₂ = a + (1 + λ₃) * α;
        l₃ = a + (1 + λ₃) / (2 + λ₃) * α;
    else
        l₂ = l₁;
        l₃ = l₁;
    end 

    ## solving discrete model
    function uniform_points_combined(L1, L2, L3, m)
        x1 = LinRange(0,        L1,        m+1)[2:end]
        x2 = LinRange(L1,       L1+L2,     m+1)[2:end]
        x3 = LinRange(L1+L2,    L1+L2+L3, m+1)[2:end]
        return vcat(collect(x1), collect(x2), collect(x3))
    end

    Domain = MCTG.DomainProperties_t(N=3, m = 10, domain_type = "1D")
    CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t()
    Embed = MCTG.CellEvent_t()
    ProlifEmbed = MCTG.CellEvent_t()
    SimTime = MCTG.SimTime_t(Tmax=maxT, δt=0.001, event_δt=0.001)
    IC = uniform_points_combined(l₁, l₂, l₃, Domain.m)
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
    # case where normalized cell length doest = 1.0 at the last node.
    if y0[end-1] == 0.0
        y0[end-1] = Cell_densities[end]
    end
    y0[end] = L0
    tspan = (0.0, maxT)
    # solve PDE
    prob = ODEProblem(MCTG.rhs!, y0, tspan, p)
    sol_cont = solve(prob, Rodas5P(), saveat=vcat([0.0:0.01:maxT]...))
    println("PDE solved")

    # Create animation 
    fig = Figure(size=(1300, 1400), title="m = $(Domain.m), N = $(Domain.N)")
    ax_springs = Axis(fig[1, 1:2], xlabel=L"x", limits=(0, maximum([L0 + 0.5, 16] ), -0.1, 0.1))
    ax_density = Axis(fig[2, 1], xlabel=L"x", ylabel=L"\rho(x,t)", limits=(0, maximum([L0 + 0.5, 16]), 0, maximum([maximum(sol_disc.Density[1]) + 0.1, 0.4])))
    ax_L_t = Axis(fig[2, 2], xlabel=L"t", ylabel=L"L(t)", limits=(0, 5, minimum([L0 - 0.5, 14]), maximum([L0 + 0.5, 16])))

    # precacluting L(t) for discrete and continuum models
    L_disc = [sol_disc.u[i][1, end] for i in eachindex(sol_disc.t)]
    L_cont = [sol_cont.u[i][end] for i in eachindex(sol_cont.t)]

    anim_name = animation_name * "special_type_$special_type.gif"
    record(fig, anim_name, 1:501; framerate=20) do frame
        t_idx = frame
        t_val = sol_cont.t[t_idx]
        t_val_title = @sprintf("t = %.2f", t_val)
        # Discrete spring model
        empty!(ax_springs)
        ax_springs.title = "N = $(Domain.N), m = $(Domain.m), k = $k, η = $η, a = $a"
        x_springs = sol_disc.u[t_idx][1,:]
        y_springs = zeros(length(x_springs))
        lines!(ax_springs, x_springs, y_springs, linewidth=3, color=:black)
        scatter!(ax_springs, x_springs, y_springs, markersize=20, marker=:circle, color=:grey)
        scatter!(ax_springs, x_springs[1:Domain.m:end], y_springs[1:Domain.m:end], markersize=20, marker=:circle, color=:red)
        text!(ax_springs, 8.0, 0.08, text = t_val_title, align = (:center, :center), fontsize=48)

        # Density evolution
        empty!(ax_density)
        disc_density = sol_disc.Density[t_idx]
        mid_disc_spring_pos = (x_springs[1:end-1] + x_springs[2:end]) ./ 2
        #ax.limits = (0, sol.u[end][end], 0, maximum(sol.u[1][1:end-1]) * 1.1)
        scatter!(ax_density, mid_disc_spring_pos, disc_density, markersize=20, marker=:circle, color=:red, label="Discrete")
        #scatter!(ax_density, mid_disc_spring_pos[1:Domain.m:end], disc_density[1:Domain.m:end], markersize=20, marker=:circle, color=:red)
        lines!(ax_density, collect(LinRange(0, 1.0, p.N)) .* sol_cont.u[t_idx][end], sol_cont.u[t_idx][1:end-1], 
            color=:blue, linewidth=3, label="Continuum")
        axislegend(ax_density)

        # Evolution of spring lengths
        empty!(ax_L_t)
        lines!(ax_L_t, sol_disc.t[1:t_idx], L_disc[1:t_idx], label="Discrete", linewidth=3, color=:red)
        lines!(ax_L_t, sol_cont.t[1:t_idx], L_cont[1:t_idx], label="Continuum", linewidth=3, color=:blue)
        if L0 > Domain.N * a
            axislegend(ax_L_t, position=:rt)
        else
            axislegend(ax_L_t, position=:rb)
        end
    end
    
    return nothing

end

k = 1; a = 5.0; η = 1; l₁ = 3.0;
generate_3_spring_FB_animation(k, a, η, l₁, "FB_3_spring_compression_anim"; special_type="None")
generate_3_spring_FB_animation(k, a, η, l₁, "FB_3_spring_compression_C1_anim"; special_type="C1")
generate_3_spring_FB_animation(k, a, η, l₁, "FB_3_spring_compression_C2_anim"; special_type="C2")
generate_3_spring_FB_animation(k, a, η, l₁, "FB_3_spring_compression_C3_anim"; special_type="C3")

k = 1; a = 5.0; η = 1; l₁ = 7.0; 
generate_3_spring_FB_animation(k, a, η, l₁, "FB_3_spring_tension_anim"; special_type="None")
generate_3_spring_FB_animation(k, a, η, l₁, "FB_3_spring_tension_C1_anim"; special_type="C1")
generate_3_spring_FB_animation(k, a, η, l₁, "FB_3_spring_tension_C2_anim"; special_type="C2")
generate_3_spring_FB_animation(k, a, η, l₁, "FB_3_spring_tension_C3_anim"; special_type="C3")

function generate_spring_FB_animation(k, a, η, L0, animation_name; maxT=10.0, special_type="None")

    Domain = MCTG.DomainProperties_t(N=10, m = 1, domain_type = "1D")
    CellMech = MCTG.CellMechProperties_t(kₛ=k, kf = 0, a = a, restoring_force="hookean")
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t()
    Embed = MCTG.CellEvent_t()
    ProlifEmbed = MCTG.CellEvent_t()
    SimTime = MCTG.SimTime_t(Tmax=maxT, δt=0.001, event_δt=0.001)
    IC = collect(LinRange(0, L0, Domain.N * Domain.m + 1)[2:end])
    NumSaveTimePoints = 1001
    sol_disc = MCTG.FreeBoundarySimulation_given_IC(IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, 1, NumSaveTimePoints);


    # Create animation 
    fig = Figure(size=(1800, 980), title="m = $(Domain.m), N = $(Domain.N)")
    ax_springs = Axis(fig[1, 1], xlabel=L"x", limits=(0, maximum([L0 + 0.5, CellMech.a * Domain.N + 1] ), -0.1, 0.1))
    ax_L_t = Axis(fig[2, 1], xlabel=L"t", ylabel=L"L(t)", limits=(0, maxT/2, minimum([L0 - 0.5, 14]), maximum([L0 + 0.5, CellMech.a * Domain.N + 1] )))

    # precacluting L(t) for discrete and continuum models
    L_disc = [sol_disc.u[i][1, end] for i in eachindex(sol_disc.t)]

    anim_name = animation_name * "special_type_$special_type.gif"
    record(fig, anim_name, 1:501; framerate=20) do frame
        t_idx = frame
        t_val = sol_disc.t[t_idx]
        t_val_title = @sprintf("t = %.2f", t_val)
        # Discrete spring model
        empty!(ax_springs)
        ax_springs.title = "N = $(Domain.N), k = $k, η = $η, a = $a"
        x_springs = sol_disc.u[t_idx][1,:]
        y_springs = zeros(length(x_springs))
        lines!(ax_springs, x_springs, y_springs, linewidth=5, color=:black)
        #scatter!(ax_springs, x_springs, y_springs, markersize=25, marker=:circle, color=:grey)
        scatter!(ax_springs, x_springs[1:Domain.m:end], y_springs[1:Domain.m:end], markersize=25, marker=:circle, color=:red)
        text!(ax_springs, 25.5, 0.08, text = t_val_title, align = (:center, :center), fontsize=48)

        # Evolution of spring lengths
        empty!(ax_L_t)
        lines!(ax_L_t, sol_disc.t[1:t_idx], L_disc[1:t_idx], label="Discrete", linewidth=5, color=:black)
    end
    
    return nothing

end

k= 4; η = 1; a = 5; L0 = 10.0; T = 50.0;       # eta*  (avoiding name collision with `eta`)
generate_spring_FB_animation(k, a, η, L0, "FB_spring_anim"; maxT=T, special_type="None")