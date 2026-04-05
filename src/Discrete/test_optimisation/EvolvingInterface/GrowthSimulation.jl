
# =====================================================================
#  Shared simulation setup helpers (eliminates code duplication)
# =====================================================================

"""Compute save-time grid with correct decimal precision."""
function _compute_save_times(SimTime, NumSaveTimePoints)
    savetimes = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    num_digits = length(string(SimTime.δt)) - 2
    return floor.(savetimes, digits=num_digits)
end

"""Create standard callback set for growth simulations."""
function _build_callbacks(SimTime, st, state)
    event_cb = PeriodicCallback(event_affect!, SimTime.δt; save_positions=(false, false))

    saved_embed_count = SavedValues(Float64, Float64)
    save_embed_count_cb = SavingCallback(store_embedded_cell_count, saved_embed_count, saveat=st)

    saved_embed_rates = SavedValues(Float64, Float64)
    save_embed_rates_cb = SavingCallback(store_embed_rates, saved_embed_rates, saveat=st)

    saved_CellMech = SavedValues(Float64, Any)
    save_CellMech_cb = SavingCallback(store_CellMech, saved_CellMech, saveat=st)

    terminate_area_cb = DiscreteCallback(area_lim_condition, terminate_affect!)
    terminate_interface_overlap_cb = DiscreteCallback(interface_overlap_condition, terminate_affect!)

    cbs = CallbackSet(save_CellMech_cb, event_cb, save_embed_count_cb,
                      save_embed_rates_cb, terminate_area_cb, terminate_interface_overlap_cb)

    return cbs, saved_embed_count, saved_embed_rates, saved_CellMech
end

"""Run the solver and collect results."""
function _solve_and_collect(prob, SimTime, st, Seed, Domain, CellMech, state,
                            saved_embed_count, saved_embed_rates, saved_CellMech, cbs, p)
    Set_Random_Seed(Seed)
    @time sol = solve(prob, Euler(), save_everystep=false, saveat=st,
                      dt=SimTime.δt, dtmax=SimTime.δt, callback=cbs, progress=true)

    embedded_cells_count = Float64[floor.(saved_embed_count.saveval)...]
    mean_embed_rates     = copy(saved_embed_rates.saveval)

    if sol.t[end] < SimTime.Tmax
        final_data = (collect(CellMech.kₛ), collect(CellMech.a),
                      collect(CellMech.kf), collect(CellMech.growth_dir))
        push!(saved_CellMech.saveval, final_data)
        if !isempty(state.embedded_cells)
            push!(embedded_cells_count,
                  floor(length(state.embedded_cells) / (Domain.m + 1)))
        end
    end

    println("Simulation Run")

    results = postSimulation(sol, p, saved_CellMech.saveval)

    if isempty(state.embedded_cells)
        embedded_matrix = Matrix{Float64}[]
    else
        embedded_matrix = convert_matrix(hcat(state.embedded_cells...), Domain.m + 1)
    end

    return results, embedded_matrix, [embedded_cells_count], [mean_embed_rates]
end


# =====================================================================
#  Main simulation entry points
# =====================================================================

"""
    GrowthSimulation(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)

Execute a 2D tissue growth simulation.
"""
function GrowthSimulation(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    M  = Int(Domain.m * Domain.N)
    st = _compute_save_times(SimTime, NumSaveTimePoints)

    state = SimulationState()
    cbs, sec, ser, scm = _build_callbacks(SimTime, st, state)

    prob, p = SetupODEproblem(M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, state)

    return _solve_and_collect(prob, SimTime, st, Seed, Domain, CellMech, state,
                              sec, ser, scm, cbs, p)
end


"""
    GrowthSimulation_with_init(...)

Execute a growth simulation with a mechanical-relaxation initialisation phase.
"""
function GrowthSimulation_with_init(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    u0 = init_problem(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)

    M  = Int(Domain.m * Domain.N)
    st = _compute_save_times(SimTime, NumSaveTimePoints)

    state = SimulationState()
    cbs, sec, ser, scm = _build_callbacks(SimTime, st, state)

    prob, p = SetupODEproblem(M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, state)
    prob = remake(prob; u0=u0)

    return _solve_and_collect(prob, SimTime, st, Seed, Domain, CellMech, state,
                              sec, ser, scm, cbs, p)
end


"""
    GrowthSimulation_given_IC(...)

Execute a growth simulation from a given initial condition.
"""
function GrowthSimulation_given_IC(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints, IC)
    M  = Int(Domain.m * Domain.N)
    st = _compute_save_times(SimTime, NumSaveTimePoints)

    state = SimulationState()
    cbs, sec, ser, scm = _build_callbacks(SimTime, st, state)

    prob, p = SetupODEproblem(M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, state)
    prob = remake(prob; u0=IC)

    return _solve_and_collect(prob, SimTime, st, Seed, Domain, CellMech, state,
                              sec, ser, scm, cbs, p)
end


"""
    init_problem(...)

Run a short mechanical-relaxation phase to generate a good initial condition.
"""
function init_problem(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    M = Int(Domain.m * Domain.N)
    CellMech_init = CellMechProperties_t(kₛ=500, growth_dir=GD_INWARD, a=20.0, kf=0)
    SimTime_init  = SimTime_t(Tmax=75, δt=0.0002, event_δt=0.0002)

    state = SimulationState()
    prob, _ = SetupODEproblem(M, Domain,
                              generate_homogeneous_population(CellMech_init, Domain.N, Domain.m),
                              SimTime_init, Prolif, Death, Embed, ProlifEmbed, state)
    @time sol = solve(prob, Euler(), save_everystep=false, dt=SimTime_init.δt, dtmax=SimTime_init.δt)
    println("Initial Problem Setup Complete")
    return sol.u[end]
end
