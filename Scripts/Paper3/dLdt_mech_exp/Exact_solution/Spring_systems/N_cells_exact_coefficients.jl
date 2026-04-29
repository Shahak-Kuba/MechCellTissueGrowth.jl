
θ  = (p, N, m) -> (2p - 1) * π / (2*N*m + 1)
λ  = (p, N, m) -> 2cos(θ(p, N, m)) - 2
vN = (p, N, m) -> cos((N*m - 0.5) * θ(p, N, m)) / cos(0.5 * θ(p, N, m))

C  = (a, p, N, m, l₀) -> (4cos.(0.5 .* θ.(p, N, m)) ./ (2N*m .+ 1)) *
                        sum(m.*cos.((j .- 0.5) .* θ.(p, N, m)) .* (l₀[j] .- a) for j in 1:N*m)

A = (a, N, l₀) -> (2 ./ N) .* (l₀ .- a)

Lₙ = (k, η, t, a, N, m, l₀) ->
    a .+ sum(C(a, p, N, m, l₀) .* vN(p, N, m) .* exp.((k .* λ(p, N, m) .* t) ./ η) for p in 1:N*m)


using CairoMakie

# figure to compare eigenvalues 
k = 1; a = 5.0; η = 1;
N = 10;
m_values = collect(1:2:50);
p_values = [1,5,10];
k = 1; a = 5.0; η = 1;
l₀ = 3.0;

#########################################
f = Figure(size=(800,600))
ax = Axis(f[1,1], aspect=1, xlabel=L"m", ylabel=L"m^{2}\lambda_{p}(m)", title = "N = $N")

for p in p_values
    asymptotic_value = -((2p - 1)^2 * π^2) / (4N^2)
    lines!(ax, m_values, asymptotic_value * ones(size(m_values)), linewidth=3, color=:black, linestyle=:dash)
    m2λ_p = zeros(length(m_values))
    for (idx, m) in enumerate(m_values)
        m2λ_p[idx] = m.^2 .* λ(p, N, m) 
    end
    scatter!(ax, m_values, m2λ_p, markersize=10, label="p = $p")
end
Legend(f[1,2], ax)
display(f)
save("exponent_asymptotic_convergence.png", f)

#########################################
N = 50;
m_values = [1,10,100];
f1 = Figure(size=(1600,1200))
ax = Axis(f1[1,1], aspect=1, xlabel=L"p/(Nm)", ylabel=L"A_{p}\mathrm{e}^{\alpha(m)}", title = L"$N = 50, \; t = 0.0$")
for m in m_values
    p_values = collect(1:2:(N.*m));
    asymptotic_value = A(a, N, fill(l₀, N*m))[p_values[1]]
    lines!(ax, p_values./ (N.*m), asymptotic_value * ones(size(p_values)), linewidth=3, color=:black, linestyle=:dash)
    A_p_exp = zeros(length(p_values))
    A_p_asymptotic_exp = zeros(length(p_values))
    for (idx, p) in enumerate(p_values)
        A_p_exp[idx] = C(a, p, N, m, fill(l₀, N*m)) * vN(p, N, m) * exp((k * m^2 * λ(p, N, m) * 0) / η)
        A_p_asymptotic_exp[idx] = A(a, N, fill(l₀, N*m))[p] * exp((k * m^2 * λ(p, N, m) * 0) / η)
    end
    println("Number of values greater than 1e-2 for m = $m: ", count(x->(abs(x)>1e-3), A_p_exp))
    scatter!(ax, p_values ./ (N.*m), A_p_exp, markersize=15, label="m = $m")
    scatter!(ax, p_values ./ (N.*m), A_p_asymptotic_exp, markersize=10, color=:red)
end
ax = Axis(f1[1,2], aspect=1, xlabel=L"p/(Nm)", ylabel=L"A_{p}\mathrm{e}^{\alpha(m)}", title = L"$N = 50, \; t = 0.1$")
for m in m_values
    p_values = collect(1:2:(N.*m));
    asymptotic_value = A(a, N, fill(l₀, N*m))[p_values[1]]
    lines!(ax, p_values./ (N.*m), asymptotic_value * ones(size(p_values)), linewidth=3, color=:black, linestyle=:dash)
    A_p_exp = zeros(length(p_values))
    A_p_asymptotic_exp = zeros(length(p_values))
    for (idx, p) in enumerate(p_values)
        A_p_exp[idx] = C(a, p, N, m, fill(l₀, N*m)) * vN(p, N, m) * exp((k * m^2 * λ(p, N, m) * 0.1) / η)
        A_p_asymptotic_exp[idx] = A(a, N, fill(l₀, N*m))[p] * exp((k * m^2 * λ(p, N, m) * 0.1) / η)
    end
    println("Number of values greater than 1e-2 for m = $m: ", count(x->(abs(x)>1e-3), A_p_exp))
    scatter!(ax, p_values ./ (N.*m), A_p_exp, markersize=15, label="m = $m")
    scatter!(ax, p_values ./ (N.*m), A_p_asymptotic_exp, markersize=10, color=:red)
end
ax = Axis(f1[2,1], aspect=1, xlabel=L"p/(Nm)", ylabel=L"A_{p}\mathrm{e}^{\alpha(m)}", title = L"$N = 50, \; t = 1.0$")
for m in m_values
    p_values = collect(1:2:(N.*m));
    asymptotic_value = A(a, N, fill(l₀, N*m))[p_values[1]]
    lines!(ax, p_values./ (N.*m), asymptotic_value * ones(size(p_values)), linewidth=3, color=:black, linestyle=:dash)
    A_p_exp = zeros(length(p_values))
    A_p_asymptotic_exp = zeros(length(p_values))
    for (idx, p) in enumerate(p_values)
        A_p_exp[idx] = C(a, p, N, m, fill(l₀, N*m)) * vN(p, N, m) * exp((k * m^2 * λ(p, N, m) * 1) / η)
        A_p_asymptotic_exp[idx] = A(a, N, fill(l₀, N*m))[p] * exp((k * m^2 * λ(p, N, m) * 1) / η)
    end
    println("Number of values greater than 1e-2 for m = $m: ", count(x->(abs(x)>1e-3), A_p_exp))
    scatter!(ax, p_values ./ (N.*m), A_p_exp, markersize=15, label="m = $m")
    scatter!(ax, p_values ./ (N.*m), A_p_asymptotic_exp, markersize=10, color=:red)
end
ax = Axis(f1[2,2], aspect=1, xlabel=L"p/(Nm)", ylabel=L"A_{p}\mathrm{e}^{\alpha(m)}", title = L"$N = 50, \; t = 10.0$")
for m in m_values
    p_values = collect(1:2:(N.*m));
    asymptotic_value = A(a, N, fill(l₀, N*m))[p_values[1]]
    lines!(ax, p_values./ (N.*m), asymptotic_value * ones(size(p_values)), linewidth=3, color=:black, linestyle=:dash)
    A_p_exp = zeros(length(p_values))
    A_p_asymptotic_exp = zeros(length(p_values))
    for (idx, p) in enumerate(p_values)
        A_p_exp[idx] = C(a, p, N, m, fill(l₀, N*m)) * vN(p, N, m) * exp((k * m^2 * λ(p, N, m) * 10) / η)
        A_p_asymptotic_exp[idx] = A(a, N, fill(l₀, N*m))[p] * exp((k * m^2 * λ(p, N, m) * 10) / η)
    end
    println("Number of values greater than 1e-2 for m = $m: ", count(x->(abs(x)>1e-3), A_p_exp))
    scatter!(ax, p_values ./ (N.*m), A_p_exp, markersize=15, label="m = $m")
    scatter!(ax, p_values ./ (N.*m), A_p_asymptotic_exp, markersize=10, color=:red)
end
axislegend(ax, position=:rb)
display(f1)
save("eigenvalues_m_values.png", f1)


#########################################
# looking at time scaling for a non-zero p coefficients
N = 50; m_vals = [100];
t_array = LinRange(0, 10, 100)
f_0 = Figure(size=(1200,1200))
∇_array = []

p = 3;
ax = Axis(f_0[1,1], aspect=1, ylabel=L"\left|A_{p}\text{exp}\left[-\frac{k^{*}m^{2}\lambda_{p}(m)t}{\eta}\right]\right|", title = "p = $p")#, limits=(-0.1, 10.1, -0.08, 0.01))
for m in m_vals
    A_p_exp = zeros(length(t_array))
    for (idx, t) in enumerate(t_array)
        A_p_exp[idx] = C(a, p, N, m, fill(l₀, N*m)) * vN(p, N, m) * exp((k * m^2 * λ(p, N, m) * t) / η)
    end
    scatter!(ax, t_array, log.(abs.(A_p_exp)), markersize=10, label="m = $m")
    ∇ = diff(log.(abs.(A_p_exp))) ./ diff(t_array)
    println(∇[1])
    push!(∇_array, ∇[1])
end

p = 10;
ax = Axis(f_0[1,2], aspect=1, title = "p = $p")#, limits=(-0.1, 10.1, -0.08, 0.01))
for m in m_vals
    A_p_exp = zeros(length(t_array))
    for (idx, t) in enumerate(t_array)
        A_p_exp[idx] = C(a, p, N, m, fill(l₀, N*m)) * vN(p, N, m) * exp((k * m^2 * λ(p, N, m) * t) / η)
    end
    scatter!(ax, t_array, log.(abs.(A_p_exp)), markersize=10, label="m = $m")
    ∇ = diff(log.(abs.(A_p_exp))) ./ diff(t_array)
    println(∇[1])
    push!(∇_array, ∇[1])
end

p = 30;
ax = Axis(f_0[2,1], aspect=1, xlabel=L"t", ylabel=L"\left|A_{p}\text{exp}\left[-\frac{k^{*}m^{2}\lambda_{p}(m)t}{\eta}\right]\right|", title = "p = $p")#, limits=(-0.1, 10.1, -0.08, 0.01))
for m in m_vals
    A_p_exp = zeros(length(t_array))
    for (idx, t) in enumerate(t_array)
        A_p_exp[idx] = C(a, p, N, m, fill(l₀, N*m)) * vN(p, N, m) * exp((k * m^2 * λ(p, N, m) * t) / η)
    end
    scatter!(ax, t_array, log.(abs.(A_p_exp)), markersize=10, label="m = $m")
    ∇ = diff(log.(abs.(A_p_exp))) ./ diff(t_array)
    println(∇[1])
    push!(∇_array, ∇[1])
end

p = 40;
ax = Axis(f_0[2,2], aspect=1, xlabel=L"t", title = "p = $p")#, limits=(-0.1, 10.1, -0.08, 0.01))
for m in m_vals
    A_p_exp = zeros(length(t_array))
    for (idx, t) in enumerate(t_array)
        A_p_exp[idx] = C(a, p, N, m, fill(l₀, N*m)) * vN(p, N, m) * exp((k * m^2 * λ(p, N, m) * t) / η)
    end
    scatter!(ax, t_array, log.(abs.(A_p_exp)), markersize=10, label="m = $m")
    ∇ = diff(log.(abs.(A_p_exp))) ./ diff(t_array)
    println(∇[1])
    push!(∇_array, ∇[1])
end
axislegend(ax, position=:rb)
display(f_0)

# looking at slope of log(|Aₚ exp[α(m)t]|) vs p for different m values
m = 1;
p = collect(1:N*m);
∇_array = []
for (idx, p_val) in enumerate(p)
    A_p_exp = zeros(length(t_array))
    for (idxx, t) in enumerate(t_array)
        A_p_exp[idxx] = C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η)
    end
    ∇ = diff(log.(abs.(A_p_exp))) ./ diff(t_array)
    println("p = $p_val, ∇ = ", ∇[1])
    push!(∇_array, ∇[1])
end

f_01 = Figure(size=(600,600))
ax = Axis(f_01[1,1], aspect=1, xlabel=L"p", ylabel=L"\nabla", title = "Slope of log(|Aₚ exp[α(m)t]|) vs p")
scatter!(ax, p, ∇_array, markersize=15)
display(f_01)


#########################################
m = 1; t = 0.1; N = 50;
p = collect(1:50);
f2 = Figure(size=(1800,600))
ax = Axis(f2[1,1], aspect=1, xlabel=L"p", ylabel=L"\left|A_{p}\text{exp}\left[-\frac{k^{*}m^{2}\lambda_{p}(m)t}{\eta}\right]\right|", title = "m = $m, t = $t")
magnitutes = zeros(length(p))
for (idx, p_val) in enumerate(p)
    magnitutes[idx] = abs(C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η))
end
barplot!(ax, p, magnitutes, color=:blue)
t = 1;
ax = Axis(f2[1,2], aspect=1, xlabel=L"p", title = "m = $m, t = $t")
magnitutes = zeros(length(p))
for (idx, p_val) in enumerate(p)
    magnitutes[idx] = abs(C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η))
end
barplot!(ax, p, magnitutes, color=:blue)
t = 10;
ax = Axis(f2[1,3], aspect=1, xlabel=L"p", title = "m = $m, t = $t")
magnitutes = zeros(length(p))
for (idx, p_val) in enumerate(p)
    magnitutes[idx] = abs(C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η))
end
barplot!(ax, p, magnitutes, color=:blue)
display(f2)
save("eigenmodes_m_$m.png", f2)


#########################################
m = 5; t = 0.1; N = 50;
p = collect(1:50);
f3 = Figure(size=(1800,600))
ax = Axis(f3[1,1], aspect=1, xlabel=L"p", ylabel=L"\left|A_{p}\text{exp}\left[-\frac{k^{*}m^{2}\lambda_{p}(m)t}{\eta}\right]\right|", title = "m = $m, t = $t")
magnitutes = zeros(length(p))
for (idx, p_val) in enumerate(p)
    magnitutes[idx] = abs(C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η))
end
barplot!(ax, p, magnitutes, color=:blue)
t = 1;
ax = Axis(f3[1,2], aspect=1, xlabel=L"p", title = "m = $m, t = $t")
magnitutes = zeros(length(p))
for (idx, p_val) in enumerate(p)
    magnitutes[idx] = abs(C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η))
end
barplot!(ax, p, magnitutes, color=:blue)
t = 10;
ax = Axis(f3[1,3], aspect=1, xlabel=L"p", title = "m = $m, t = $t")
magnitutes = zeros(length(p))
for (idx, p_val) in enumerate(p)
    magnitutes[idx] = abs(C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η))
end
barplot!(ax, p, magnitutes, color=:blue)
display(f3)
save("eigenmodes_m_$m.png", f3)


#########################################
m = 100; t = 0.1; N = 50;
p = collect(1:50);
f4 = Figure(size=(1800,600))
ax = Axis(f4[1,1], aspect=1, xlabel=L"p", ylabel=L"\left|A_{p}\text{exp}\left[-\frac{k^{*}m^{2}\lambda_{p}(m)t}{\eta}\right]\right|", title = "m = $m, t = $t")
magnitutes = zeros(length(p))
for (idx, p_val) in enumerate(p)
    magnitutes[idx] = abs(C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η))
end
barplot!(ax, p, magnitutes, color=:blue)
t = 1;
ax = Axis(f4[1,2], aspect=1, xlabel=L"p", title = "m = $m, t = $t")
magnitutes = zeros(length(p))
for (idx, p_val) in enumerate(p)
    magnitutes[idx] = abs(C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η))
end
barplot!(ax, p, magnitutes, color=:blue)
t = 10;
ax = Axis(f4[1,3], aspect=1, xlabel=L"p", title = "m = $m, t = $t")
magnitutes = zeros(length(p))
for (idx, p_val) in enumerate(p)
    magnitutes[idx] = abs(C(a, p_val, N, m, fill(l₀, N*m)) * vN(p_val, N, m) * exp((k * m^2 * λ(p_val, N, m) * t) / η))
end
barplot!(ax, p, magnitutes, color=:blue)
display(f4)
save("eigenmodes_m_$m.png", f4)


#########################################
f5 = Figure(size=(650,600))
ax = Axis(f5[1,1], aspect=1, xlabel=L"m", ylabel=L"A_{p}", title = "N = $N")
m_vals = collect(1:2:50);
p = [1, 10, 50];
A_p_asymptotic = A(a, N, fill(l₀, N*m))[p[1]]
lines!(ax, m_vals, A_p_asymptotic .* ones(size(m_vals)), linewidth=3, color=:black, linestyle=:dash, label = "Asym value")
for p_val in p
    magnitutes_A = zeros(length(m_vals))
    for (idx, m_val) in enumerate(m_vals)
        magnitutes_A[idx] = C(a, p_val, N, m_val, fill(l₀, N*m_val)) * vN(p_val, N, m_val)
    end
    scatter!(ax, m_vals, magnitutes_A, markersize=10, label="p = $p_val")
end
axislegend(ax, position=:rt)
display(f5)
save("A_p_convergence.png", f5)


#########################################
f6 = Figure(size=(1800,600))
m_vals = collect(1:2:50);
p = [1, 10, 50];
ϵ = [0.1, 0.01, 0.001];
for (kk, ϵ_val) in enumerate(ϵ)
    if kk == 1
        ax = Axis(f6[1,kk], aspect=1, xlabel=L"m", ylabel=L"A_{p}", title = "ϵ = $ϵ_val")
    else
        ax = Axis(f6[1,kk], aspect=1, xlabel=L"m", title = "ϵ = $ϵ_val")
    end
    A_p_asymptotic = A(a, N, fill(l₀, N*m))[p[1]]
    lines!(ax, m_vals, A_p_asymptotic .* ones(size(m_vals)), linewidth=3, color=:black, linestyle=:dash, label = "Asym value")
    for p_val in p
        MCTG.Set_Random_Seed(1)
        magnitutes_A = zeros(length(m_vals))
        for (idx, m_val) in enumerate(m_vals)
            l_array = zeros(N*m_val)
            for ii in eachindex(l_array)
                l_array[ii] = l₀ + ϵ_val*rand()
            end
            #println(l_array)
            magnitutes_A[idx] = C(a, p_val, N, m_val, l_array) * vN(p_val, N, m_val)
        end
        scatter!(ax, m_vals, magnitutes_A, markersize=10, label="p = $p_val")
    end
    if kk == 3
        axislegend(ax, position=:rt)
    end
end
display(f6)
save("A_p_convergence_random.png", f6)


#########################################
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
