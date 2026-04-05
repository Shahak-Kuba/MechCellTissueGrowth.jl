import MechCellTissueGrowth as MCTG
using JLD2

output_file = "/Users/shahakmac/Desktop/PhD/Code/MechCellTissueGrowth.jl/outData/Feb_2026/BoneSimData_ks5_LB.jld2"
@load output_file all_solutions all_embedded_cell_pos embedded_count_iteration_results Ω_iteration_results

ξ_all = zeros(size(embedded_count_iteration_results[1],1)-1, size(embedded_count_iteration_results,1));
A_all = zeros(size(embedded_count_iteration_results[1],1)-1, size(embedded_count_iteration_results,1));
N_embed_all = zeros(size(embedded_count_iteration_results[1],1)-1, size(embedded_count_iteration_results,1));
for iteration in 1: size(all_solutions,1)
    ξ,A,N_embed = MCTG.WallDensity(embedded_count_iteration_results[iteration], Ω_iteration_results[iteration])
    for ii in eachindex(ξ)
        if ξ[ii] < 0.0
            ξ[ii] = 0.0
        end
    end
    if length(ξ_all[:, iteration]) == length(ξ)
        ξ_all[:, iteration] .= ξ
    else
        n = min(length(ξ_all[:, iteration]), length(ξ))
        ξ_all[1:n, iteration] .= ξ[1:n]
        ξ_all[n+1:end, iteration] .= NaN
    end
    if length(A_all[:, iteration]) == length(A)
        A_all[:, iteration] .= A
    else
        n = min(length(A_all[:, iteration]), length(A))
        A_all[1:n, iteration] .= A[1:n]
        A_all[n+1:end, iteration] .= NaN
    end
    if length(N_embed_all[:, iteration]) == length(N_embed)
        N_embed_all[:, iteration] .= N_embed
    else
        n = min(length(N_embed_all[:, iteration]), length(N_embed))
        N_embed_all[1:n, iteration] .= N_embed[1:n]
        N_embed_all[n+1:end, iteration] .= NaN
    end
end
idx = 1
sol = all_solutions[idx];
embedded_pos = all_embedded_cell_pos[idx];
fig7_1, fig7_2 = MCTG.PlotWallEmbedDensities(sol, embedded_pos, ξ_all, 0.000625)
display(fig7_1)

# Plot 8: Bone formation WTh_asymmetry_ratio
results = MCTG.analyse_wallThickness_density(all_solutions, all_embedded_cell_pos);
sort!(results, :WThRatio);
println(results)
fig8_1, fig8_2 = MCTG.PlotWThAsymmetry(results, all_solutions, all_embedded_cell_pos);

min_Vn, mean_Vn, max_Vn = MCTG.AnalyseBoneFormationVelocity(all_solutions)

#save figures
save(joinpath(pwd(), "Scripts/Paper2/RawFigures/fig7_1.png"), fig7_1)
save(joinpath(pwd(), "Scripts/Paper2/RawFigures/fig7_2.png"), fig7_2)
save(joinpath(pwd(), "Scripts/Paper2/RawFigures/fig8_1.png"), fig8_1)
save(joinpath(pwd(), "Scripts/Paper2/RawFigures/fig8_2.png"), fig8_2)
