"""
    sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,prolif,death,embed,α,β,γ,event_δt,seed, NumSaveTimePoints)

Execute a series of 2D mechanical relaxation simulations.

This function sets up and runs 2D simulations for different stiffness coefficients and boundary types. It iterates over arrays of stiffness coefficients and boundary types, sets up the corresponding ODE problem for each case, solves it using a specified numerical method with periodic callbacks, and collects the results.

# Simulation Parameters
- `N`: Number of cells in the simulation.
- `m`: Number of springs per cell.
- `R₀`: Radius or characteristic length of the initial shape.
- `D`: Array of diffuision coefficients used to calculate cell stiffness.
- `l₀`: Resting length of the spring per cell.
- `kf`: Tissue production rate per cell.
- `η`: Viscosity or damping coefficient per cell.
- `growth_dir`: Direction of tissue growth ('inward' or 'outward').
- `Tmax`: Total simulation time (in days).
- `δt`: Time step for the numerical integration.
- `btypes`: Types of boundary conditions (e.g., 'circle', 'triangle').
- `dist_type`: Distribution type for node placement (e.g., 'Linear', 'sigmoid').
- `prolif`, `death`, `embed`: Boolean flags indicating cell behaviors.
- `α`, `β`, `γ`: Parameters for cell behaviors.
- `event_δt`: Time interval for periodic callback events.
- `seed`: Seed number for reproducability with stochastic simulations.

# Returns
A vector of vectors of `SimResults_t` objects. Each inner vector represents the simulation results for different boundary types under a specific stiffness coefficient.

# Example
```julia
all_results = sim2D(N,m,R₀,D,l₀,kf,η,growth_dir,Tmax,δt,btypes,dist_type,
                        prolif, death, embed, α, β, γ, event_δt, seed, NumSaveTimePoints);
```
"""

function GrowthSimulation(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    M = Int(Domain.m * Domain.N) # total number of springs along the interface
    savetimes = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    # calculating how many digits are in SimTime.δt 
    num_digits = (length(string(SimTime.δt)) - 2) # -2 to remove the 0. from the string and -1 to remove the last digit
    # ensure that the save times are of the form x.xxxx
    st = floor.(savetimes, digits=num_digits)

    # for cell embedment 
    global embedded_cells = []
    global cell_embedment_times = []
    embedded_cells_count = []
    mean_embed_rates = []

    # cell prolif, death, embedment event callback
    event_cb = PeriodicCallback(event_affect!, SimTime.δt; save_positions = (false, false))
    # saving callback 
    saved_embed_count = SavedValues(Float64, Float64)
    save_embed_count_cb = SavingCallback(store_embedded_cell_count, saved_embed_count, saveat = st)
    saved_embed_rates = SavedValues(Float64, Float64)
    save_embed_rates_cb = SavingCallback(store_embed_rates, saved_embed_rates, saveat = st)
    saved_CellMech = SavedValues(Float64 , Vector{Vector{Any}})
    save_CellMech_cb = SavingCallback(store_CellMech, saved_CellMech, saveat = st)
    # Terminate simulation callback (based on area coverage)
    terminate_area_cb = DiscreteCallback(area_lim_condition, terminate_affect!)
    terminate_interface_overlap_cb = DiscreteCallback(interface_overlap_condition, terminate_affect!)
    # generating a set of callbacks
    cbs = CallbackSet(save_CellMech_cb, event_cb, save_embed_count_cb, save_embed_rates_cb, terminate_area_cb, terminate_interface_overlap_cb)

    prob, p = SetupODEproblem(M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    Set_Random_Seed(Seed)
    @time sol = solve(prob, Euler(), save_everystep = false, saveat = st, dt = SimTime.δt, dtmax = SimTime.δt, callback = cbs, progress = true)
    push!(embedded_cells_count, floor.(saved_embed_count.saveval))
    push!(mean_embed_rates, saved_embed_rates.saveval)
    if sol.t[end] < SimTime.Tmax
        final_CellMech_data = []
        push!(final_CellMech_data,CellMech.kₛ)
        push!(final_CellMech_data,CellMech.a)
        push!(final_CellMech_data,CellMech.kf)
        push!(final_CellMech_data,CellMech.growth_dir)
        push!(saved_CellMech.saveval, final_CellMech_data)
        push!(embedded_cells_count[1], floor.(size(hcat(embedded_cells...),2)/(Domain.m+1)))
    end
    println("Simulation Run")
    return postSimulation(sol, p, saved_CellMech.saveval), convert_matrix(hcat(embedded_cells...), Domain.m + 1), embedded_cells_count, mean_embed_rates
end

function init_problem(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    M = Int(Domain.m * Domain.N) # total number of springs along the interface
    CellMech_init = CellMechProperties_t(kₛ=500, growth_dir="inward", a=20.0, kf = 0)
    SimTime_init = SimTime_t(Tmax=75, δt = 0.0002, event_δt = 0.0002)
    prob, p = SetupODEproblem(M, Domain, CellMech_init, SimTime_init, Prolif, Death, Embed, ProlifEmbed)
    @time sol = solve(prob, Euler(), save_everystep = false, dt = SimTime_init.δt, dtmax = SimTime_init.δt)
    println("Initial Problem Setup Complete")
    return sol.u[end]
end

function GrowthSimulation_with_init(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints)
    u0 = init_problem(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    
    # actual tissue growth simulation
    M = Int(Domain.m * Domain.N)
    savetimes = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    # calculating how many digits are in SimTime.δt 
    num_digits = (length(string(SimTime.δt)) - 2) # -2 to remove the 0. from the string and -1 to remove the last digit
    # ensure that the save times are of the form x.xxxx
    st = floor.(savetimes, digits=num_digits)

    # for cell embedment 
    global embedded_cells = []
    global cell_embedment_times = []
    embedded_cells_count = []
    mean_embed_rates = []

    # cell prolif, death, embedment event callback
    event_cb = PeriodicCallback(event_affect!, SimTime.δt; save_positions = (false, false))
    # saving callback 
    saved_embed_count = SavedValues(Float64, Float64)
    save_embed_count_cb = SavingCallback(store_embedded_cell_count, saved_embed_count, saveat = st)
    saved_embed_rates = SavedValues(Float64, Float64)
    save_embed_rates_cb = SavingCallback(store_embed_rates, saved_embed_rates, saveat = st)
    saved_CellMech = SavedValues(Float64 , Vector{Vector{Any}})
    save_CellMech_cb = SavingCallback(store_CellMech, saved_CellMech, saveat = st)
    # Terminate simulation callback (based on area coverage)
    terminate_area_cb = DiscreteCallback(area_lim_condition, terminate_affect!)
    terminate_interface_overlap_cb = DiscreteCallback(interface_overlap_condition, terminate_affect!)
    # generating a set of callbacks
    cbs = CallbackSet(save_CellMech_cb, event_cb, save_embed_count_cb, save_embed_rates_cb, terminate_area_cb, terminate_interface_overlap_cb)

    prob, p = SetupODEproblem(M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    prob.u0 = u0
    Set_Random_Seed(Seed)
    @time sol = solve(prob, Euler(), save_everystep = false, saveat = st, dt = SimTime.δt, dtmax = SimTime.δt, callback = cbs, progress = true)
    push!(embedded_cells_count, floor.(saved_embed_count.saveval))
    push!(mean_embed_rates, saved_embed_rates.saveval)
    if sol.t[end] < SimTime.Tmax
        final_CellMech_data = []
        push!(final_CellMech_data,CellMech.kₛ)
        push!(final_CellMech_data,CellMech.a)
        push!(final_CellMech_data,CellMech.kf)
        push!(final_CellMech_data,CellMech.growth_dir)
        push!(saved_CellMech.saveval, final_CellMech_data)
        push!(embedded_cells_count[1], floor.(size(hcat(embedded_cells...),2)/(Domain.m+1)))
    end
    println("Simulation Run")
    return postSimulation(sol, p, saved_CellMech.saveval), convert_matrix(hcat(embedded_cells...), Domain.m + 1), embedded_cells_count, mean_embed_rates
end


function GrowthSimulation_given_IC(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, NumSaveTimePoints, IC)
    
    # actual tissue growth simulation
    M = Int(Domain.m * Domain.N)
    savetimes = LinRange(0, SimTime.Tmax, NumSaveTimePoints)
    # calculating how many digits are in SimTime.δt 
    num_digits = (length(string(SimTime.δt)) - 2) # -2 to remove the 0. from the string and -1 to remove the last digit
    # ensure that the save times are of the form x.xxxx
    st = floor.(savetimes, digits=num_digits)

    # for cell embedment 
    global embedded_cells = []
    global cell_embedment_times = []
    embedded_cells_count = []
    mean_embed_rates = []

    # cell prolif, death, embedment event callback
    event_cb = PeriodicCallback(event_affect!, SimTime.δt; save_positions = (false, false))
    # saving callback 
    saved_embed_count = SavedValues(Float64, Float64)
    save_embed_count_cb = SavingCallback(store_embedded_cell_count, saved_embed_count, saveat = st)
    saved_embed_rates = SavedValues(Float64, Float64)
    save_embed_rates_cb = SavingCallback(store_embed_rates, saved_embed_rates, saveat = st)
    saved_CellMech = SavedValues(Float64 , Vector{Vector{Any}})
    save_CellMech_cb = SavingCallback(store_CellMech, saved_CellMech, saveat = st)
    # Terminate simulation callback (based on area coverage)
    terminate_area_cb = DiscreteCallback(area_lim_condition, terminate_affect!)
    terminate_interface_overlap_cb = DiscreteCallback(interface_overlap_condition, terminate_affect!)
    # generating a set of callbacks
    cbs = CallbackSet(save_CellMech_cb, event_cb, save_embed_count_cb, save_embed_rates_cb, terminate_area_cb, terminate_interface_overlap_cb)

    prob, p = SetupODEproblem(M, Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
    prob.u0 = IC
    Set_Random_Seed(Seed)
    @time sol = solve(prob, Euler(), save_everystep = false, saveat = st, dt = SimTime.δt, dtmax = SimTime.δt, callback = cbs, progress = true)
    push!(embedded_cells_count, floor.(saved_embed_count.saveval))
    push!(mean_embed_rates, saved_embed_rates.saveval)
    if sol.t[end] < SimTime.Tmax
        final_CellMech_data = []
        push!(final_CellMech_data,CellMech.kₛ)
        push!(final_CellMech_data,CellMech.a)
        push!(final_CellMech_data,CellMech.kf)
        push!(final_CellMech_data,CellMech.growth_dir)
        push!(saved_CellMech.saveval, final_CellMech_data)
        push!(embedded_cells_count[1], floor.(size(hcat(embedded_cells...),2)/(Domain.m+1)))
    end
    println("Simulation Run")
    return postSimulation(sol, p, saved_CellMech.saveval), convert_matrix(hcat(embedded_cells...), Domain.m + 1), embedded_cells_count, mean_embed_rates
end