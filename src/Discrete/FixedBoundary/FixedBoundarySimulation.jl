function FixedBoundarySimulation(FB_IC, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    M = Int(Domain.m * Domain.N) # total number of springs along the interface
    savetimes = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    # calculating how many digits are in SimTime.δt 
    num_digits = (length(string(SimTime.δt)) - 2) # -2 to remove the 0. from the string and -1 to remove the last digit
    # ensure that the save times are of the form x.xxxx
    st = floor.(savetimes, digits=num_digits)

    # for cell embedment 
    #global embedded_cells = []
    #global cell_embedment_times = []
    #embedded_cells_count = []
    #mean_embed_rates = []

    # cell prolif, death, embedment event callback
    event_cb = PeriodicCallback(event_affect!, SimTime.event_δt; save_positions = (false, false))
    # saving callback 
    #saved_embed_count = SavedValues(Float64, Float64)
    #save_embed_count_cb = SavingCallback(store_embedded_cell_count, saved_embed_count, saveat = st)
    #saved_embed_rates = SavedValues(Float64, Float64)
    #save_embed_rates_cb = SavingCallback(store_embed_rates, saved_embed_rates, saveat = st)
    saved_CellMech = SavedValues(Float64 , Vector{Vector{Any}})
    save_CellMech_cb = SavingCallback(store_CellMech, saved_CellMech, saveat = st)
    # Terminate simulation callback (based on area coverage)
    #terminate_area_cb = DiscreteCallback(area_lim_condition, terminate_affect!)
    #terminate_interface_overlap_cb = DiscreteCallback(interface_overlap_condition, terminate_affect!)
    # generating a set of callbacks
    cbs = CallbackSet(save_CellMech_cb, event_cb)

    HomCellMech = generate_homogeneous_population_FB(CellMech, Domain.N, Domain.m);
    #println("ks = " * string(HomCellMech.kₛ[1]) * ", η = " * string(HomCellMech.η) * ", a = " * string(HomCellMech.a[1]))
    prob, p = SetupFixedBoundaryODEproblem(FB_IC, M, Domain, HomCellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    
    Set_Random_Seed(Seed)
   
    @time sol = solve(prob, Euler(), save_everystep = false, saveat = st, dt = SimTime.δt, dtmax = SimTime.δt, callback = cbs, progress = true)
   
    #push!(embedded_cells_count, floor.(saved_embed_count.saveval))
    #push!(mean_embed_rates, saved_embed_rates.saveval)
    if sol.t[end] < SimTime.Tmax
        final_CellMech_data = []
        push!(final_CellMech_data,CellMech.kₛ)
        push!(final_CellMech_data,CellMech.a)
        push!(final_CellMech_data,CellMech.kf)
        push!(final_CellMech_data,CellMech.growth_dir)
        push!(saved_CellMech.saveval, final_CellMech_data)
        #push!(embedded_cells_count[1], floor.(size(hcat(embedded_cells...),2)/(Domain.m+1)))
    end

    println("Simulation Run")
    return postSimulation(sol, p, saved_CellMech.saveval)
end


function FixedBoundarySimulation_given_IC(IC::Vector{Float64}, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    M = Int(Domain.m * Domain.N) # total number of springs along the interface
    savetimes = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    # calculating how many digits are in SimTime.δt 
    num_digits = (length(string(SimTime.δt)) - 2) # -2 to remove the 0. from the string and -1 to remove the last digit
    # ensure that the save times are of the form x.xxxx
    st = floor.(savetimes, digits=num_digits)

    # for cell embedment 
    #global embedded_cells = []
    #global cell_embedment_times = []
    #embedded_cells_count = []
    #mean_embed_rates = []

    # cell prolif, death, embedment event callback
    event_cb = PeriodicCallback(event_affect!, SimTime.event_δt; save_positions = (false, false))
    # saving callback 
    saved_CellMech = SavedValues(Float64 , Vector{Vector{Any}})
    save_CellMech_cb = SavingCallback(store_CellMech, saved_CellMech, saveat = st)
    # generating a set of callbacks
    cbs = CallbackSet(save_CellMech_cb, event_cb)

    HomCellMech = generate_homogeneous_population_FB(CellMech, Domain.N, Domain.m);
    println("k = $(HomCellMech.kₛ), η = $(HomCellMech.η)")
    #println("ks = " * string(HomCellMech.kₛ[1]) * ", η = " * string(HomCellMech.η) * ", a = " * string(HomCellMech.a[1]))
    prob, p = SetupFixedBoundaryODEproblem_given_IC(IC, M, Domain, HomCellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    
    Set_Random_Seed(Seed)
   
    @time sol = solve(prob, Euler(), save_everystep = false, saveat = st, dt = SimTime.δt, dtmax = SimTime.δt, callback = cbs, progress = true)
   
    if sol.t[end] < SimTime.Tmax
        final_CellMech_data = []
        push!(final_CellMech_data,CellMech.kₛ)
        push!(final_CellMech_data,CellMech.a)
        push!(final_CellMech_data,CellMech.kf)
        push!(final_CellMech_data,CellMech.growth_dir)
        push!(saved_CellMech.saveval, final_CellMech_data)
        #push!(embedded_cells_count[1], floor.(size(hcat(embedded_cells...),2)/(Domain.m+1)))
    end

    println("Simulation Run")
    return postSimulation(sol, p, saved_CellMech.saveval)
end