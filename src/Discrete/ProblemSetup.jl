function generate_homogeneous_population(C1::CellMechProperties_t, N, m)
    M = N * m
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a = ElasticArray{Float64}(zeros(M))
    growth_dir = ElasticArray{String}([])

    for i in 1:M
        ks[i] = C1.kₛ * m
        kf[i] = C1.kf / m
        a[i] = C1.a / m
        push!(growth_dir, C1.growth_dir)
    end

    return HeterogeneousCellMechProperties_t(ks, C1.η / m, C1.restoring_force, kf, growth_dir, a)
end

function generate_heterogenous_population(C1::CellMechProperties_t, C2::CellMechProperties_t, N, m, C1_count)
    M = N * m
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a = ElasticArray{Float64}(zeros(M))
    growth_dir = ElasticArray{String}([])

    for i in 1:M
        if i <= C1_count * m
            ks[i] = C1.kₛ * m
            kf[i] = C1.kf / m
            a[i] = C1.a / m
            push!(growth_dir, C1.growth_dir)
        else
            ks[i] = C2.kₛ * m
            kf[i] = C2.kf / m
            a[i] = C2.a / m
            push!(growth_dir, C2.growth_dir)
        end
    end

    return HeterogeneousCellMechProperties_t(ks, C1.η / m, C1.restoring_force, kf, growth_dir, a)
end


function generate_heterogenous_population(C1::CellMechProperties_t, C2::CellMechProperties_t, C3::CellMechProperties_t, N, m, C1_count, C2_count)
    M = N * m
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a = ElasticArray{Float64}(zeros(M))
    growth_dir = ElasticArray{String}([])

    for i in 1:M
        if i <= C1_count * m
            ks[i] = C1.kₛ * m
            kf[i] = C1.kf / m
            a[i] = C1.a / m
            push!(growth_dir, C1.growth_dir)
        elseif C1_count * m < i <= (C1_count + C2_count) * m
            ks[i] = C2.kₛ * m
            kf[i] = C2.kf / m
            a[i] = C2.a / m
            push!(growth_dir, C2.growth_dir)
        else
            ks[i] = C3.kₛ * m
            kf[i] = C3.kf / m
            a[i] = C3.a / m
            push!(growth_dir, C3.growth_dir)
        end
    end

    return HeterogeneousCellMechProperties_t(ks, C1.η / m, C1.restoring_force, kf, growth_dir, a)
end

function generate_heterogenous_population(CellTypes::Vector{CellMechProperties_t}, N, m, Cell_counts::Vector{Int})
    # make sure that cell counts add up to n
    if sum(Cell_counts) != N
        throw(ArgumentError("Cell counts must add up to N"))
    end
    M = N * m
    ks = ElasticArray{Float64}(zeros(M))
    kf = ElasticArray{Float64}(zeros(M))
    a = ElasticArray{Float64}(zeros(M))
    growth_dir = ElasticArray{String}([])

    # Need to generalise this code!
    for ii in 1:length(CellTypes)
        if ii == 1
            jj = 1
        else
            jj = sum(Cell_counts[1:ii-1] .* m) + 1
        end
        println("ii = ", ii," jj = ",jj)
        while jj <= sum(Cell_counts[1:ii] .* m)
            ks[jj] = CellTypes[ii].kₛ * m
            kf[jj] = CellTypes[ii].kf / m
            a[jj] = CellTypes[ii].a / m
            push!(growth_dir, CellTypes[ii].growth_dir)
            jj += 1
        end
    end

    return HeterogeneousCellMechProperties_t(ks, CellTypes[1].η / m, CellTypes[1].restoring_force, kf, growth_dir, a)
end

"""
    SetupODEproblem(btype, M, m, R₀, kₛ, η, kf, l₀, δt, Tmax, growth_dir, prolif, death, embed, α, β, γ, dist_type)

Set up and configure a 1D or 2D ODE problem for mechanical relaxation simulations in tissue growth.

This function initializes the conditions and parameters for a 2D ODE problem based on the specified boundary type, physical parameters, and cell behaviors. It then constructs an ODEProblem object, ready for solving with DifferentialEquations.jl.

# Arguments
- `btype`: Boundary type (e.g., 'circle', 'triangle').
- `M`: Total number of springs along the interface.
- `m`: Number of springs per cell.
- `R₀`: Initial radius or characteristic length of the shape.
- `kₛ`: Spring stiffness coefficient.
- `η`: Viscosity or damping coefficient.
- `kf`: Tissue production rate per cell.
- `l₀`: Resting length of the spring per cell.
- `δt`: Time step for the numerical integration.
- `Tmax`: Total simulation time.
- `growth_dir`: Direction of tissue growth ('inward' or 'outward').
- `prolif`: Boolean flag for cell proliferation.
- `death`: Boolean flag for cell death.
- `embed`: Boolean flag for cell embedding.
- `α, β, γ`: Parameters for the cell behaviors.
- `dist_type`: Distribution type for node placement.

# Returns
- `ODEProblem`: An ODE problem instance set up with the specified parameters and initial conditions.
- `p`: A tuple containing the parameters used in setting up the ODE problem.
"""

function SetupODEproblem(M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    #aₘ = CellMech.a ./ Domain.m
    #kₘ = CellMech.kₛ .* Domain.m
    #ηₘ = CellMech.η ./ Domain.m
    #kfₘ = CellMech.kf ./ Domain.m
    u0 = u0SetUp(Domain.btype,Domain.R₀,M,Domain.dist_type,Domain.domain_type)
    #p = (Domain.m,kₛ,η,kf,l₀,SimTime.δt,growth_dir,domain_type,btype,restoring_force,Prolif,Death,Embed,SimTime)
    p = (Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    tspan = (0.0, SimTime.Tmax)
    return ODEProblem(Growth_ODE!,u0,tspan,p), p
end
