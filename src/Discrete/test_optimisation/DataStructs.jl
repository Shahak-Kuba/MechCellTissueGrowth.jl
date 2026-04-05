
# ============================= Enums =============================
"""Restoring force types for cell springs."""
@enum ForceType Hookean Nonlinear Nonlinear2

"""Growth direction: stored as Int8 so the sign can be used directly in normal vector computation."""
@enum GrowthDir::Int8 GD_INWARD=-1 GD_OUTWARD=1

"""Event trigger modes for cell behaviour callbacks."""
@enum EventTrigger ET_Constant ET_Periodic

"""Event probability function types."""
@enum EventFunc EF_Constant EF_Constant2 EF_Length

# ============================= Results =============================
"""
    SimResults_t

A data structure representing the results of a simulation.

# Fields
- `btype`: A `String` indicating the type of simulation or boundary condition.
- `t`: Time steps at which the simulation results are recorded.
- `u`: State of the system at each time step.
- `∑F`: Sum of forces at each time step.
- `Density`: Density of the system at each time step.
- `Vₙ`: Normal velocity at each time step.
- `Ω`: Void area within the pore at each time step.
- `σ`: Stress of the springs/cells at each time step.
- `Κ`: Approximate curvature of the moving boundary at each time step.
- `CellCount`: Amount of active cells at a given time `t`.
"""
struct SimResults_t
    btype::String
    t::Vector{Float64}
    u::Vector{Matrix{Float64}}
    ∑F::Vector{Vector{Float64}}
    Density::Vector{Vector{Float64}}
    Vₙ::Vector{Vector{Float64}}
    Ω::Vector{Float64}
    σ::Vector{Vector{Float64}}
    Κ::Vector{Vector{Float64}}
    CellCount::Vector{Int64}
end

# ============================= Cell Event =============================
"""
    CellEvent_t

Encapsulates the properties related to cell events in a simulation.

### Fields
- `flag::Bool`: Whether the cell event is active. (Default: `false`)
- `rate::Float64`: The rate at which the event occurs. (Default: `0.0`)
- `event_func::EventFunc`: The function that defines the event's behaviour. (Default: `EF_Constant`)
"""
@kwdef struct CellEvent_t
    flag::Bool = false
    rate::Float64 = 0.0
    event_func::EventFunc = EF_Constant
end

# ============================= Cell Mechanical Properties =============================
"""
    CellMechProperties_t

Holds the mechanical properties of a cell (homogeneous — all cells share these values).

### Fields
- `kₛ::Float64`: Spring stiffness. (Default: `25.0`)
- `η::Float64`: Damping coefficient. (Default: `1.0`)
- `restoring_force::ForceType`: Restoring force type. (Default: `Nonlinear`)
- `kf::Float64`: Tissue production rate per cell. (Default: `87.842`)
- `growth_dir::GrowthDir`: Direction of cell growth. (Default: `GD_INWARD`)
- `D::Float64`: Diffusion coefficient. (Default: `0.0`)
- `a::Float64`: Rest length parameter. (Default: `10.0`)
"""
@kwdef struct CellMechProperties_t
    kₛ::Float64 = 25.0
    η::Float64 = 1.0
    restoring_force::ForceType = Nonlinear
    kf::Float64 = 87.842
    growth_dir::GrowthDir = GD_INWARD
    D::Float64 = 0.0
    a::Float64 = 10.0
end

"""
    HeterogeneousCellMechProperties_t

Per-spring mechanical properties for heterogeneous cell populations.
`growth_dir` is stored as `Int8` (-1 = inward, +1 = outward) for direct use as a sign in the normal vector computation.
"""
@kwdef mutable struct HeterogeneousCellMechProperties_t
    kₛ::ElasticArray{Float64}
    η::Float64 = 1.0
    restoring_force::ForceType = Nonlinear
    kf::ElasticArray{Float64}
    growth_dir::ElasticArray{Int8}
    a::ElasticArray{Float64}
end

# ============================= Domain Properties =============================
"""
    DomainProperties_t

Defines the properties of the simulation domain.

### Fields
- `N::Int64`: Total number of cells. (Default: `100`)
- `m::Int64`: Number of springs per cell. (Default: `6`)
- `R₀::Float64`: Shape radius in μm. (Default: `105.0`)
- `domain_type::String`: `"1D"` or `"2D"`. (Default: `"2D"`)
- `btype::String`: Shape of the domain. (Default: `"circle"`)
- `dist_type::String`: Distribution type for cells. (Default: `"Linear"`)
"""
@kwdef struct DomainProperties_t
    N::Int64 = 100
    m::Int64 = 6
    R₀::Float64 = 105.0
    domain_type::String = "2D"
    btype::String = "circle"
    dist_type::String = "Linear"
    Ω_0::Float64 = π * R₀^2
    dir_to_EXP_Image::String = ""
end

# ============================= Simulation Time =============================
"""
    SimTime_t

Time-related properties for the simulation.

### Fields
- `Tmax::Float64`: Maximum simulation time in days. (Default: `20.0`)
- `δt::Float64`: Time step. (Default: `0.01`)
- `event_trigger::EventTrigger`: Method for triggering events. (Default: `ET_Constant`)
- `event_δt::Float64`: Time interval for triggering events. (Default: `0.01`)
- `periodic_δt::Float64`: Period for periodic events. (Default: `4.0`)
- `event_length::Float64`: Duration of periodic event window. (Default: `0.0`)
"""
@kwdef struct SimTime_t
    Tmax::Float64 = 20.0
    δt::Float64 = 0.01
    event_trigger::EventTrigger = ET_Constant
    event_δt::Float64 = 0.01
    periodic_δt::Float64 = 4.0
    event_length::Float64 = 0.0
end

# ============================= ODE Cache =============================
"""
    ODECache

Pre-allocated buffers for the ODE right-hand side, eliminating per-timestep allocations.
All vectors are length N (number of spring nodes).
"""
mutable struct ODECache
    N::Int
    l::Vector{Float64}       # spring lengths
    ρ::Vector{Float64}       # spring densities (1/l)
    dvx::Vector{Float64}     # directional vector x-component
    dvy::Vector{Float64}     # directional vector y-component
    τx::Vector{Float64}      # tangent vector x-component
    τy::Vector{Float64}      # tangent vector y-component
    nx::Vector{Float64}      # normal vector x-component
    ny::Vector{Float64}      # normal vector y-component
    F_tang::Vector{Float64}  # tangential force projection
end

function ODECache(N::Int)
    ODECache(N,
        zeros(N), zeros(N),
        zeros(N), zeros(N),
        zeros(N), zeros(N),
        zeros(N), zeros(N),
        zeros(N)
    )
end

function resize_cache!(cache::ODECache, N::Int)
    cache.N = N
    for field in (:l, :ρ, :dvx, :dvy, :τx, :τy, :nx, :ny, :F_tang)
        resize!(getfield(cache, field), N)
    end
    nothing
end

# ============================= Simulation State =============================
"""
    SimulationState

Replaces global variables for tracking embedded cells during simulation.
"""
mutable struct SimulationState
    embedded_cells::Vector{Vector{Float64}}
    cell_embedment_times::Vector{Float64}
end

SimulationState() = SimulationState(Vector{Float64}[], Float64[])

# ============================= Free Boundary IC =============================
@kwdef struct FB_IC_t
    q0::Function = x -> 2.0
    q0_der::Function = x -> 0.0
    L0::Float64 = 10.0
end

struct ContinuumSolution_t
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
end
