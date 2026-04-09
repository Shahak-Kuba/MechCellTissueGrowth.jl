using LinearAlgebra
using Printf
using CairoMakie

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

function eigenvalue_near_zero(Matrix_func::Function, M::Int, k::Int=1)
    A = Matrix_func(M)
    vals = eigvals(Symmetric(A))   # ascending order
    return vals[end-k+1]
end

function AnalyseEigenvalsOfSystem(Matrix_func::Function, plt_title::String)
    # Matrix sizes
    Ms = [10, 20, 50, 100, 250, 500]

    # Smallest-in-magnitude eigenvalue
    lam1 = [eigenvalue_near_zero(Matrix_func, M, 1) for M in Ms]

    x = log.(Float64.(Ms))
    y = log.(abs.(lam1))
    X = hcat(ones(length(x)), x)
    coef = X \ y
    intercept, slope = coef

    # Fitted values
    fit_vals = exp(intercept) .* Float64.(Ms).^slope

    # Rescaled values: should approach π^2
    #scaled_vals = [(M+1)^2 * abs(λ) for (M, λ) in zip(Ms, lam1)]

    # -----------------------------
    # Figure 1: log-log scaling plot
    # -----------------------------
    fig1 = Figure(size = (700, 700))
    ax1 = Axis(
        fig1[1, 1],
        aspect = 1,
        title = plt_title,
        xlabel = L"\log (M)",
        ylabel = L"\log (|λ_{1}(M)|)",
        xscale = log10,
        yscale = log10
    )

    lines!(ax1, Ms, abs.(lam1), label = L"\log(|λ_{1}(M)|)", linewidth = 3)
    scatter!(ax1, Ms, abs.(lam1), marker='*', markersize=40)
    lines!(ax1, Ms, fit_vals, linestyle = :dash, linewidth = 3, label = @sprintf("fit slope = %.2f", slope))
    axislegend(ax1, position = :rt)

    return fig1

end

fixed_boundary_fig = AnalyseEigenvalsOfSystem(tri_matrix_dense_fixed, "Fixed boundaries");
display(fixed_boundary_fig)

free_boundary_fig = AnalyseEigenvalsOfSystem(tri_matrix_dense_free, "Free boundary");
display(free_boundary_fig)

# explicit form for eigenvalues
λₚ_fixed = (p,M) -> -2 .+ 2 .* cos.((π * p)/ M)
λₚ_free = (p,M) -> -2 .+ 2 .* cos.((π * p)/ (M+1)) + (2/(M+1)).*(sin.((π * p)/ (M+1))).^2

function compare_numeric_analytic_eigenvals(Matrix_func::Function, eigenval_func::Function, plt_title::String, p::Int64=1)
    # Matrix sizes
    Ms = [10, 20, 50, 100, 250, 500]

    lam_numeric = [eigenvalue_near_zero(Matrix_func, M, p) for M in Ms]
    lam_analytic = [eigenval_func(p, M) for M in Ms]
    println(lam_analytic - lam_numeric)

    fig1 = Figure(size = (1200, 500))
    ax1 = Axis(
        fig1[1, 1],
        aspect = 1,
        title = plt_title,
        xlabel = L"M",
        ylabel = L"λ_{1}(M)"
    )
    lines!(ax1, Ms, lam_numeric, linewidth = 3, label="Numerical");
    scatter!(ax1, Ms, lam_numeric, marker='*',  markersize=40)
    lines!(ax1, Ms, lam_analytic, linewidth = 3, label="Analytical");
    scatter!(ax1, Ms, lam_analytic, marker='*',  markersize=40)
    axislegend(ax1, position = :rb)

    ax2 = Axis(
        fig1[1, 2],
        aspect = 1,
        xlabel = L"M",
        ylabel = L"λ_{1}^{(\text{A})}(M) - λ_{1}^{(\text{N})}(M)",
        title = "Difference"
    )
    lines!(ax2, Ms, lam_analytic - lam_numeric, linewidth = 3);
    scatter!(ax2, Ms, lam_analytic -  lam_numeric, marker='*',  markersize=40)

    return fig1
end

fixed_boundary_fig_2 = compare_numeric_analytic_eigenvals(tri_matrix_dense_fixed, λₚ_fixed, "Fixed boundaries");
display(fixed_boundary_fig_2)

free_boundary_fig_2 = compare_numeric_analytic_eigenvals(tri_matrix_dense_free, λₚ_free, "Free boundary");
display(free_boundary_fig_2)

save("fixed_boundary_loglog.png", fixed_boundary_fig)
save("free_boundary_loglog.png", free_boundary_fig)
save("fixed_boundary_Numerical_vs_Analytical.png", fixed_boundary_fig_2)
save("free_boundary_Numerical_vs_Analytical.png", free_boundary_fig_2)
