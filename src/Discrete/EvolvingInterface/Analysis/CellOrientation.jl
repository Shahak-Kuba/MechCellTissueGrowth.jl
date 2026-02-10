function cell_orientation_plot(embedded_cells, t)
    cell_orientation = Float64[]
    for cell in embedded_cells
        cell_left_edge = cell[:,1]
        cell_right_edge = cell[:,end]
        cell_length = KubaPhD.δ(cell_left_edge', cell_right_edge')[1]
        θ_x = acos((cell_right_edge[1] - cell_left_edge[1]) / cell_length)
        push!(cell_orientation, θ_x)
    end

    unique_cell_orientation = unique(round.(cell_orientation, digits=2))
    freq = Int64[]
    for val in unique_cell_orientation
        push!(freq, count(i -> (i == val), round.(cell_orientation, digits=2)))
    end

    θ_care_range = 0:π/12:π
    freq_by_θ_care_range = []
    for ii in 2:size(θ_care_range, 1)
        v = 0
        for jj in eachindex(unique_cell_orientation)
            if (θ_care_range[ii-1] <= unique_cell_orientation[jj]) && (unique_cell_orientation[jj] < θ_care_range[ii])
                v += freq[jj]
            end
        end
        push!(freq_by_θ_care_range, v)
    end

    θ_rep = repeat(θ_care_range, inner=2)
    freq_rep = repeat(freq_by_θ_care_range, inner=2)

    θ = vcat(θ_rep, θ_rep[3:end] .+ π)[2:end]
    R = vcat(freq_rep ./ 2, freq_rep ./ 2, freq_rep[1] / 2)

    text_t = floor(t)

    Plots.plot!(θ, R, proj=:polar, line=:step, linewidth=4,color=:blue, title="t = $text_t", legend=false)
end

function time_orientation_plot(all_embedded_cells, embedded_cell_count, idx, t)
    P = Plots.plot()
    cell_orientation_plot(reverse(all_embedded_cells[1:Int64(embedded_cell_count[1][idx])]), t[idx])
    return P
end


#function time_orientation_plot_all(all_embedded_cells, embedded_cell_count, t)
#    P = Plots.plot()
#    for ii in size(embedded_cell_count[1],1):-1:2
#        cell_orientation_plot(reverse(all_embedded_cells[1:Int64(embedded_cell_count[1][ii])]), t)
#    end
#    return P
#end