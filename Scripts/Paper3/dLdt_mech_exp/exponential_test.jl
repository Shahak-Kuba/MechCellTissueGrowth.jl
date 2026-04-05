using CairoMakie

A = 2.0; B = 5.0; λ₁ = 5.0; λ₂ = 1.0;
m = 2;

f = t -> (A/m) * exp(-λ₁*m*t) + (B) * exp(-λ₂*t) 
t_vals = 0.0:0.001:5.0

F = Figure(size=(800, 400))
ax1 = Axis(F[1,1], xlabel=L"t", ylabel=L"f(t)", title="m = $(m)")
lines!(ax1, t_vals, f.(t_vals), linewidth=3, color=:black)
ax2 = Axis(F[1,2], xlabel=L"t", ylabel=L"\log (f(t))", title="m = $(m)")
lines!(ax2, t_vals[1:100], log.(f.(t_vals[1:100])), linewidth=3, color=:black)
display(F)

t1 = 0.0; t2 = 0.01;
∇ = (f(t2) - f(t1)) / (t2 - t1)
println("λ₁ ≈ ", ∇/m)