import MechCellTissueGrowth as MCTG
using JLD2
using SMTPClient
using Dates

BatchSize = 1
dir_to_exp_image = pwd() * "/EXP_IMAGES/FM40-2-S2-SGH-Outline.png"
E_law = "Constant"


function run_bone_simulations(BatchSize, dir_to_exp_image, E_law)
    # Domain Properties
    Domain = MCTG.DomainProperties_t(N=45, R₀=115.0, btype="Exp_Image", m = 5, dir_to_EXP_Image=dir_to_exp_image)

    # Cell Mechanics
    CellMech = MCTG.CellMechProperties_t(kₛ=5, growth_dir="outward", kf = 40)

    # Cell Behaviours
    Prolif = MCTG.CellEvent_t()
    Death = MCTG.CellEvent_t(false, 0.0001, "Constant")
    Embed = MCTG.CellEvent_t(true, 0.000625*CellMech.kf, E_law)
    ProlifEmbed = MCTG.CellEvent_t()
    # Simulation time parameters
    SimTime = MCTG.SimTime_t(Tmax=70, δt=0.001, event_δt=0.001)

    function run_simulation_with_init(CellMech, Domain, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, IC)
        HomCellMech = MCTG.generate_homogeneous_population(CellMech, Domain.N, Domain.m);
        num_t_save = 5#round(Int64, SimTime.Tmax / 2) 
        sol, embed_pos, embed_count, embed_rates = MCTG.GrowthSimulation_given_IC(Domain, HomCellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, num_t_save, IC);
        return sol, embed_pos, embed_count

    end

    embedded_count_iteration_results = Vector{Int64}[]
    Ω_iteration_results = Vector{Float64}[]
    t = Vector{Float64}[]
    all_solutions = MCTG.SimResults_t[]
    all_embedded_cell_pos = Vector{Matrix{Float64}}[]

    prob_u0 = []
    MAX_ITERS = 2 * BatchSize

    iter = 1
    successful_simulations = 1

    while iter < MAX_ITERS && successful_simulations < BatchSize + 1

        if iter == 1
            prob_u0 = MCTG.init_problem(Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed)
        end
        #HomCellMech = MCTG.generate_homogeneous_population(CellMech, Domain.N, Domain.m);
        sol, embedded_cells, embed_cell_count = run_simulation_with_init(CellMech, Domain, SimTime, Prolif, Death, Embed, ProlifEmbed, iter, prob_u0)

        push!(all_solutions, sol)
        push!(all_embedded_cell_pos, embedded_cells)
        push!(embedded_count_iteration_results, convert(Vector{Int64}, embed_cell_count[1]))
        push!(Ω_iteration_results, sol.Ω[1] .- sol.Ω)
        if successful_simulations == 1
            push!(t, sol.t)
        end
        println("Simulation $successful_simulations / $BatchSize")
        successful_simulations += 1

        iter += 1
    end

    return all_solutions, all_embedded_cell_pos, embedded_count_iteration_results, Ω_iteration_results, t
end

all_solutions, all_embedded_cell_pos, embedded_count_iteration_results, Ω_iteration_results, t = run_bone_simulations(BatchSize, dir_to_exp_image, E_law)

idx = 1
sol = all_solutions[idx];
embedded_cells = all_embedded_cell_pos[idx];
MCTG.PlotOsteon_Simulation(sol, embedded_cells)