import MechCellTissueGrowth as MCTG
using JLD2
using SMTPClient
using Dates

BatchSize = 3
# Domain Properties
Domain = MCTG.DomainProperties_t(N=45, R₀=60.0, btype="PerturbedCircle", m = 5)

# Cell Mechanics
CellMech = MCTG.CellMechProperties_t(kₛ=5, growth_dir="inward", kf = 35)

# Cell Behaviours
Prolif = MCTG.CellEvent_t()
Death = MCTG.CellEvent_t(true, 0.0001, "Constant")
Embed = MCTG.CellEvent_t(true, 0.000625*CellMech.kf, "Constant")
ProlifEmbed = MCTG.CellEvent_t()
# Simulation time parameters
#SimTime = MCTG.SimTime_t(Tmax=12, event_trigger="", event_δt=8.0, event_length=2.0)
SimTime = MCTG.SimTime_t(Tmax=56, δt=0.001, event_δt=0.001)

function run_simulation_with_init(CellMech, Domain, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, IC)
    HomCellMech = MCTG.generate_homogeneous_population(CellMech, Domain.N, Domain.m);
    num_t_save = round(Int64, SimTime.Tmax / 7) 
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
    HomCellMech = MCTG.generate_homogeneous_population(CellMech, Domain.N, Domain.m);
    sol, embedded_cells, embed_cell_count = run_simulation_with_init(CellMech, Domain, SimTime, Prolif, Death, Embed, ProlifEmbed, iter, prob_u0)

    if sol.t[end] == SimTime.Tmax
        push!(all_solutions, sol)
        push!(all_embedded_cell_pos, embedded_cells)
        push!(embedded_count_iteration_results, convert(Vector{Int64}, embed_cell_count[1]))
        push!(Ω_iteration_results, sol.Ω[1] .- sol.Ω)
        if successful_simulations == 1
            push!(t, sol.t)
        end
        println("Simulation $successful_simulations / $BatchSize")
        successful_simulations += 1
    else
        println("Simulation $iter failed to reach Tmax. Final time: ", sol.t[end])
    end

    iter += 1
end

MCTG.send_email_alert("shahakjuliacode@gmail.com", "jwdz aiva nuxa nnpc", "Bone Batch Simulations Complete", "s.kuba@qut.edu.au");

# saving data
current_dir = pwd()
date_time = Dates.now()
output_dir = current_dir * "/outData/" * string(date_time)
if !isdir(output_dir)
    mkdir(output_dir)
end
filename = "/BoneBatchSimData.jld2"
output_file = output_dir * filename
@save output_file all_solutions all_embedded_cell_pos embedded_count_iteration_results Ω_iteration_results



