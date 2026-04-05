function FreeBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    M          = Int(Domain.m * Domain.N)
    savetimes  = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    num_digits = length(string(SimTime.δt)) - 2
    st         = floor.(savetimes, digits=num_digits)

    # SimulationState replaces global variables — consistent with GrowthSimulation
    state = SimulationState()

    # Cell event callback — event_affect! expects 9-element p, which SetupFBODEproblem now provides
    event_cb = PeriodicCallback(event_affect!, SimTime.event_δt; save_positions=(false, false))

    saved_CellMech    = SavedValues(Float64, Any)
    save_CellMech_cb  = SavingCallback(store_CellMech, saved_CellMech, saveat=st)

    #terminate_area_cb             = DiscreteCallback(area_lim_condition, terminate_affect!)
    #terminate_interface_overlap_cb = DiscreteCallback(interface_overlap_condition, terminate_affect!)

    cbs = CallbackSet(save_CellMech_cb, event_cb)

    HomCellMech = generate_homogeneous_population_FB(CellMech, Domain.N, Domain.m)
    prob, p     = SetupFBODEproblem(FB_IC, M, Domain, HomCellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, state)

    Set_Random_Seed(Seed)
    @time sol = solve(prob, Euler(), save_everystep=false, saveat=st,
                      dt=SimTime.δt, dtmax=SimTime.δt, callback=cbs, progress=true)

    if sol.t[end] < SimTime.Tmax
        final_data = (collect(HomCellMech.kₛ), collect(HomCellMech.a),
                      collect(HomCellMech.kf), collect(HomCellMech.growth_dir))
        push!(saved_CellMech.saveval, final_data)
    end

    println("Simulation Run")
    return postSimulation(sol, p, saved_CellMech.saveval)
end


function FreeBoundarySimulation_given_IC(IC::Vector{Float64}, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    M          = Int(Domain.m * Domain.N)
    savetimes  = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    num_digits = length(string(SimTime.δt)) - 2
    st         = floor.(savetimes, digits=num_digits)

    state = SimulationState()

    event_cb = PeriodicCallback(event_affect!, SimTime.event_δt; save_positions=(false, false))

    saved_CellMech    = SavedValues(Float64, Any)
    save_CellMech_cb  = SavingCallback(store_CellMech, saved_CellMech, saveat=st)

    #terminate_area_cb             = DiscreteCallback(area_lim_condition, terminate_affect!)
    #terminate_interface_overlap_cb = DiscreteCallback(interface_overlap_condition, terminate_affect!)

    cbs = CallbackSet(save_CellMech_cb, event_cb)

    HomCellMech = generate_homogeneous_population_FB(CellMech, Domain.N, Domain.m)
    println("k = $(HomCellMech.kₛ), η = $(HomCellMech.η)")
    prob, p = SetupFBODEproblem_given_IC(IC, M, Domain, HomCellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, state)

    Set_Random_Seed(Seed)
    @time sol = solve(prob, Euler(), save_everystep=false, saveat=st,
                      dt=SimTime.δt, dtmax=SimTime.δt, callback=cbs, progress=true)

    if sol.t[end] < SimTime.Tmax
        final_data = (collect(HomCellMech.kₛ), collect(HomCellMech.a),
                      collect(HomCellMech.kf), collect(HomCellMech.growth_dir))
        push!(saved_CellMech.saveval, final_data)
    end

    println("Simulation Run")
    return postSimulation(sol, p, saved_CellMech.saveval)
end
