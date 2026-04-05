"""
    generate_homogeneous_population(C1, N, m)

Create per-spring mechanical properties for a homogeneous cell population.
Returns a `HeterogeneousCellMechProperties_t` where every spring shares the same values.
"""
function generate_homogeneous_population(C1::CellMechProperties_t, N, m)
    M = N * m
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a  = ElasticArray{Float64}(zeros(M))
    gd_val = Int8(C1.growth_dir)   # -1 or +1
    growth_dir = ElasticArray{Int8}(fill(gd_val, M))

    for i in 1:M
        ks[i] = C1.kₛ * m
        kf[i] = C1.kf / m
        a[i]  = C1.a / m
    end

    return HeterogeneousCellMechProperties_t(ks, C1.η / m, C1.restoring_force, kf, growth_dir, a)
end

function generate_homogeneous_population_FB(C1::CellMechProperties_t, N, m)
    M = (N * m) + 1
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a  = ElasticArray{Float64}(zeros(M))
    gd_val = Int8(C1.growth_dir)
    growth_dir = ElasticArray{Int8}(fill(gd_val, M))

    for i in 1:M
        ks[i] = C1.kₛ * m
        kf[i] = C1.kf / m
        a[i]  = C1.a / m
    end

    return HeterogeneousCellMechProperties_t(ks, C1.η / m, C1.restoring_force, kf, growth_dir, a)
end

function generate_heterogenous_population(C1::CellMechProperties_t, C2::CellMechProperties_t, N, m, C1_count)
    M = N * m
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a  = ElasticArray{Float64}(zeros(M))
    growth_dir = ElasticArray{Int8}(zeros(Int8, M))

    for i in 1:M
        if i <= C1_count * m
            ks[i] = C1.kₛ * m
            kf[i] = C1.kf / m
            a[i]  = C1.a / m
            growth_dir[i] = Int8(C1.growth_dir)
        else
            ks[i] = C2.kₛ * m
            kf[i] = C2.kf / m
            a[i]  = C2.a / m
            growth_dir[i] = Int8(C2.growth_dir)
        end
    end

    return HeterogeneousCellMechProperties_t(ks, C1.η / m, C1.restoring_force, kf, growth_dir, a)
end

function generate_heterogenous_population(C1::CellMechProperties_t, C2::CellMechProperties_t, C3::CellMechProperties_t, N, m, C1_count, C2_count)
    M = N * m
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a  = ElasticArray{Float64}(zeros(M))
    growth_dir = ElasticArray{Int8}(zeros(Int8, M))

    for i in 1:M
        if i <= C1_count * m
            ks[i] = C1.kₛ * m; kf[i] = C1.kf / m; a[i] = C1.a / m
            growth_dir[i] = Int8(C1.growth_dir)
        elseif C1_count * m < i <= (C1_count + C2_count) * m
            ks[i] = C2.kₛ * m; kf[i] = C2.kf / m; a[i] = C2.a / m
            growth_dir[i] = Int8(C2.growth_dir)
        else
            ks[i] = C3.kₛ * m; kf[i] = C3.kf / m; a[i] = C3.a / m
            growth_dir[i] = Int8(C3.growth_dir)
        end
    end

    return HeterogeneousCellMechProperties_t(ks, C1.η / m, C1.restoring_force, kf, growth_dir, a)
end

function generate_heterogenous_population(CellTypes::Vector{CellMechProperties_t}, N, m, Cell_counts::Vector{Int})
    sum(Cell_counts) != N && throw(ArgumentError("Cell counts must add up to N"))

    M = N * m
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a  = ElasticArray{Float64}(zeros(M))
    growth_dir = ElasticArray{Int8}(zeros(Int8, M))

    for ii in 1:length(CellTypes)
        jj = ii == 1 ? 1 : sum(Cell_counts[1:ii-1] .* m) + 1
        jj_end = sum(Cell_counts[1:ii] .* m)
        while jj <= jj_end
            ks[jj] = CellTypes[ii].kₛ * m
            kf[jj] = CellTypes[ii].kf / m
            a[jj]  = CellTypes[ii].a / m
            growth_dir[jj] = Int8(CellTypes[ii].growth_dir)
            jj += 1
        end
    end

    return HeterogeneousCellMechProperties_t(ks, CellTypes[1].η / m, CellTypes[1].restoring_force, kf, growth_dir, a)
end


"""
    SetupODEproblem(M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, state)

Set up and configure the ODE problem for 2D tissue growth simulation.
Creates the `ODECache` and includes both cache and `state` in the parameter tuple.
"""
function SetupODEproblem(M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, state)
    u0 = u0SetUp(Domain.btype, Domain.R₀, M, Domain.dist_type, Domain.domain_type; dir_to_img=Domain.dir_to_EXP_Image)
    cache = ODECache(size(u0, 2))
    p = (Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, cache, state)
    tspan = (0.0, SimTime.Tmax)
    return ODEProblem(Growth_ODE!, u0, tspan, p), p
end


# =====================================================================
#  Free Boundary Problem Setup
# =====================================================================

function generate_discrete_IC_from_density_profile(Func, Func_Derivative, N, m, L0, eps=0.1)
    M = N * m
    integral, _ = quadgk(x -> Func(x), 0, L0)
    lambda = M / integral

    function x_interior!(F, x)
        for i in 1:length(F)
            if i == 1
                F[i] = ((1 / (x[i] - 0)) + (1 / (x[i] - x[i+1]))) / ((0 - x[i+1]) / 2) - lambda * Func_Derivative(x[i])
            elseif i == length(F)
                F[i] = ((1 / (x[i] - x[i-1])) + (1 / (x[i] - L0))) / ((x[i-1] - L0) / 2) - lambda * Func_Derivative(x[i])
            else
                F[i] = ((1 / (x[i] - x[i-1])) + (1 / (x[i] - x[i+1]))) / ((x[i-1] - x[i+1]) / 2) - lambda * Func_Derivative(x[i])
            end
        end
    end

    x_guess = collect(range(eps, stop=L0 - eps, length=M - 1))
    sol = nlsolve(x_interior!, x_guess)
    return vcat(0.0, sol.zero, L0)
end

function generate_discrete_IC_from_density_profile(Func, Func_Derivative, M, L0, eps=0.1)
    integral, _ = quadgk(x -> Func(x), 0, L0)
    lambda = M / integral

    function x_interior!(F, x)
        for i in 1:length(F)
            if i == 1
                F[i] = ((1 / (x[i] - 0)) + (1 / (x[i] - x[i+1]))) / ((0 - x[i+1]) / 2) - lambda * Func_Derivative(x[i])
            elseif i == length(F)
                F[i] = ((1 / (x[i] - x[i-1])) + (1 / (x[i] - L0))) / ((x[i-1] - L0) / 2) - lambda * Func_Derivative(x[i])
            else
                F[i] = ((1 / (x[i] - x[i-1])) + (1 / (x[i] - x[i+1]))) / ((x[i-1] - x[i+1]) / 2) - lambda * Func_Derivative(x[i])
            end
        end
    end

    x_guess = collect(range(eps, stop=L0 - eps, length=M - 1))
    sol = nlsolve(x_interior!, x_guess)
    return vcat(0.0, sol.zero, L0)
end

function SetupFBODEproblem(FB_IC, M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, state)
    q0     = FB_IC.q0
    q0_der = FB_IC.q0_der
    L0     = FB_IC.L0

    x0 = generate_discrete_IC_from_density_profile(q0, q0_der, M, L0, 0.001)
    u0 = ElasticMatrix([x0'; zeros(1, length(x0))])
    cache = ODECache(size(u0, 2))
    p = (Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, cache, state)
    tspan = (0.0, SimTime.Tmax)
    return ODEProblem(FB_ODE!, u0, tspan, p), p
end

function SetupFBODEproblem_given_IC(IC::Vector{Float64}, M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, state)
    u0 = ElasticMatrix(zeros(2, M + 1))
    u0[1, 2:end] = IC
    cache = ODECache(size(u0, 2))
    p = (Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, cache, state)
    tspan = (0.0, SimTime.Tmax)
    return ODEProblem(FB_ODE!, u0, tspan, p), p
end
