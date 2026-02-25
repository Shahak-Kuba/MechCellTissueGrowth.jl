import MechCellTissueGrowth as MCTG
using OrdinaryDiffEq

# Parameters
k = 5.0; η = 1.0; α = k/η; a = 5.0; L0 = 10.0; N = 1001;

# Force and Diffusivity functions (Hookean spring model)
Ffunc = ρ -> k * (1/ρ - a)
Dfunc = ρ -> α / ρ^2 

# Proliferation and apoptosis functions
Pfunc = ρ -> 0.0
Afunc = ρ -> 0.01

# Force and Diffusivity functions (Nonlinear spring model)
#Ffunc = ρ -> k * (1/a - ρ)
#Dfunc = ρ -> α 

p = MCTG.FBParams(α=α, η=η, k=k, a=a, L0=L0, N=N, rBC=:free, lBC=:fixed, F=Ffunc, D=Dfunc, P=Pfunc, A=Afunc)

q0_func = z -> 2.0 #* cos(z*π/2) + 1.0 #2.0


y0 = MCTG.make_initial_condition_FB(p.N; U0fun = q0_func, L0=p.L0)

tspan = (0.0, 40.0)

# solve PDE
prob = ODEProblem(MCTG.rhs!, y0, tspan, p)
sol = solve(prob, Rodas5P(), saveat=vcat([0.0:0.0001:40.0]...))

function generate_simulation_figures(sol, p, density_profile_img_name, dLdt_img_name)

    Fig = Figure(size=(650,500))
    if p.rBC == :fixed
        title_str = "Fixed boundary"
    elseif p.rBC == :free
        title_str = "Free boundary"
    end
    ax = Axis(Fig[1,1], xlabel="x", ylabel="q(x,t)", title=title_str)
    Nlines = 15
    idxs = round.(Int, range(1, length(sol.t), length = min(Nlines, length(sol.t))))
    for j in idxs
        x = range(0.0, 1.0, length=p.N) * sol.u[j][end];
        lines!(ax, x, sol.u[j][1:p.N], color=sol.t[j], colorrange = (sol.t[1], sol.t[end]), linewidth=3)
    end
    Colorbar(Fig[1,2], colorrange = (sol.t[1], sol.t[end]), label = "time", width = 15)
    save("$density_profile_img_name.png", Fig)

    Fig = Figure(size=(650,500))
    ax2 = Axis(Fig[1,1], xlabel="t", ylabel="L(t)", title="L(t) over time")
    lines!(ax2, sol.t, [sol.u[j][end] for j in 1:length(sol.t)], color=:blue)
    save("$dLdt_img_name.png", Fig)
end

generate_simulation_figures(sol, p, "qxt_nonlinear_discretisation", "Lxt_nonlinear_discretisation")
