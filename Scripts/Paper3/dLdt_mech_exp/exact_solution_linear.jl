using LinearAlgebra
# making the matrix

function tri_matrix_dense_fixed(N::Int)
    A = zeros(Float64, N, N)
    for i in 1:N
        A[i,i] = -2
        if i < N
            A[i, i+1] = 1
            A[i+1, i] = 1
        end
    end
    return A
end

function tri_matrix_dense_free(N::Int)
    A = zeros(Float64, N, N)
    for i in 1:N
        A[i,i] = -2
        if i < N
            A[i, i+1] = 1
            A[i+1, i] = 1
        end
    end
    A[N, N] = -1 # free boundary condition modification
    return A
end



k_s = 5.0; η_s = 10.0;
N = 5;
m_values = collect(1:100);

min_eigenvalues = Float64[];
max_eigenvalues = Float64[];

for m in m_values
    M = N * m;
    A = tri_matrix_dense_fixed(M);
    eigenvalues = eigvals(A);
    push!(min_eigenvalues, minimum(eigenvalues));
    push!(max_eigenvalues, maximum(eigenvalues));
end

for (i, m) in enumerate(m_values)
    M = N * m
    println("Min eigenvalue for M = $(N*m): ", min_eigenvalues[i], " Maximum eigenvalue for M = $(N*m): ", max_eigenvalues[i])
end

using CairoMakie
F = Figure(size=(900,950))
ax = Axis(F[1,1], aspect=1, xlabel=L"$M$", ylabel=L"\lambda_{\text{min}}", title="Eigenvalues of the System Matrix", xlabelsize=24, ylabelsize=24, xticklabelsize=20, yticklabelsize=20);
lines!(ax, log.(N .* m_values), min_eigenvalues, color = :blue, linewidth=4, label="Minimum eigenvalues")
#lines!(ax, N .* m_values, max_eigenvalues, color = :red, linewidth=4, label="Maximum eigenvalues")
axislegend(ax, position = :rt)
display(F)

"""
using Symbolics
M = 4;
@variables x[1:M](t) x0[1:M] t k η a
#k_c = 1.0; η_c = 1.0; a_c = 2.0;
A = tri_matrix_dense_free(M)
DiagA = eigen(A)
λ = DiagA.values
S = DiagA.vectors
#S = [1.0 -2.0 - λ[1]; 1.0 -2.0 - λ[2]]
S_inv = inv(S)
E = Matrix{Num}(undef, M, M)
for i in 1:M
    for j in 1:M
        E[i,j] = 0.0
    end
    E[i,i] = exp((k/η * λ[i] * t))
end
y_vec = zeros(Num,M); y_vec[end] = a;
z0 =  x0 + inv(A) * y_vec
x = S * E * S_inv * z0
"""


using LinearAlgebra
using Printf
using CairoMakie

# Build the M x M tridiagonal matrix
function tridiag_matrix(M::Int)
    A = zeros(Float64, M, M)
    for i in 1:M
        A[i, i] = -2.0
        if i > 1
            A[i, i-1] = 1.0
        end
        if i < M
            A[i, i+1] = 1.0
        end
    end
    return A
end

# k = 1 gives the eigenvalue closest to zero
# k = 2 gives the next one, etc.
function eigenvalue_near_zero(M::Int, k::Int=1)
    A = tridiag_matrix(M)
    vals = eigvals(Symmetric(A))   # ascending order
    return vals[end-k+1]
end

# Matrix sizes
Ms = [10, 20, 40, 80, 160, 320]

# Smallest-in-magnitude eigenvalue
lam1 = [eigenvalue_near_zero(M, 1) for M in Ms]

println(" M        lambda_1(M)              (M+1)^2 * |lambda_1(M)|")
println("----------------------------------------------------------------")
for (M, λ) in zip(Ms, lam1)
    @printf("%4d   % .10e         %.10f\n", M, λ, (M+1)^2 * abs(λ))
end

# Log-log fit: log(|λ|) = a + b log(M)
x = log.(Float64.(Ms))
y = log.(abs.(lam1))
X = hcat(ones(length(x)), x)
coef = X \ y
intercept, slope = coef

println()
@printf("Estimated log-log slope = %.6f\n", slope)

# Ratio test: should tend to 1/4 for quadratic scaling
println()
println("Ratio test: |lambda_1(2M)| / |lambda_1(M)|")
println("-------------------------------------------")
for i in 1:length(Ms)-1
    ratio = abs(lam1[i+1]) / abs(lam1[i])
    @printf("M = %3d -> %3d   ratio = %.6f\n", Ms[i], Ms[i+1], ratio)
end

# Fitted values
fit_vals = exp(intercept) .* Float64.(Ms).^slope

# Rescaled values: should approach π^2
scaled_vals = [(M+1)^2 * abs(λ) for (M, λ) in zip(Ms, lam1)]

# -----------------------------
# Figure 1: log-log scaling plot
# -----------------------------
fig1 = Figure(size = (700, 500))
ax1 = Axis(
    fig1[1, 1],
    title = "Smallest eigenvalue",
    xlabel = L"\log (M)",
    ylabel = L"\log (|λ_{1}(M)|)",
    xscale = log10,
    yscale = log10
)

lines!(ax1, Ms, abs.(lam1), label = L"\log|λ₁(M)|", linewidth = 3,)
scatter!(ax1, Ms, abs.(lam1), markersize=15)
lines!(ax1, Ms, fit_vals, linestyle = :dash, linewidth = 3, label = @sprintf("fit slope = %.2f", slope))
axislegend(ax1, position = :rt)

display(fig1)
save("eigenvalue_loglog.png", fig1)

