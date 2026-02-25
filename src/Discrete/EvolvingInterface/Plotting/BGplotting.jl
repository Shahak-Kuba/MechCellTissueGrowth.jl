function greyscale_colors(N::Integer)
    N ≥ 1 || throw(ArgumentError("N must be ≥ 1"))
    vals = range(0.3, 0.95; length = N)
    return [RGB(v, v, v) for v in vals]
end

function PlotOsteon_Simulation(sol, embedded_pos)
    fig = Figure(size=(800,800))
    ga = fig[1,1] = GridLayout()
    ax = CairoMakie.Axis(ga[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", xlabelsize = 22, xlabelvisible = false, ylabelvisible= false, ylabelsize = 22, aspect = 1, xticklabelsize = 18, yticklabelsize = 18, xticklabelsvisible=false, yticklabelsvisible=false )
    colors = reverse(KubaPhD.greyscale_colors(size(sol.t, 1)))
    
    for ii in axes(sol.t, 1)
        if ii == 1 || ii == size(sol.t, 1)
            poly!(ax, KubaPhD.convert_coordinates_to_tuples(sol.u[ii]), color=colors[ii], strokecolor = :black, strokewidth = 2, joinstyle=:round,  stroke_depth_shift = +1e-3)
        else
            poly!(ax, KubaPhD.convert_coordinates_to_tuples(sol.u[ii]), color=colors[ii], strokecolor = :grey98, strokewidth = 1, joinstyle=:round,  stroke_depth_shift = +1e-3)
        end
    end
    poly!(ax, KubaPhD.convert_coordinates_to_tuples(sol.u[end]), color=colors[1], strokecolor = :black, strokewidth = 2, joinstyle=:round)

    for cell in embedded_pos
        cell = KubaPhD.equidistant_normal_points_matrix(cell, 1.0, 15.0)
        #CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=3, color=KubaPhD.ColorsCBfriendly[:Black], joinstyle=:round)
        CairoMakie.poly!(ax, cell[1,:], cell[2,:], color=KubaPhD.ColorsCBfriendly[:Black], strokecolor=:black, strokewidth=2, joinstyle=:round, stroke_depth_shift = +1e-3)
        CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=3, color=KubaPhD.ColorsCBfriendly[:Black], joinstyle=:round)
    end
    return fig
end

function PlotOsteonSim!(ax, sol, embedded_pos)
    colors = reverse(greyscale_colors(size(sol.t, 1)))
    
    for ii in axes(sol.t, 1)
        if ii == 1 || ii == size(sol.t, 1)
            poly!(ax, convert_coordinates_to_tuples(sol.u[ii]), color=colors[ii], strokecolor = :black, strokewidth = 2, joinstyle=:round,  stroke_depth_shift = +1e-3)
        else
            poly!(ax, convert_coordinates_to_tuples(sol.u[ii]), color=colors[ii], strokecolor = :grey98, strokewidth = 1, joinstyle=:round,  stroke_depth_shift = +1e-3)
        end
    end
    poly!(ax, convert_coordinates_to_tuples(sol.u[end]), color=colors[1], strokecolor = :black, strokewidth = 2, joinstyle=:round)

    for (i,cell) in enumerate(embedded_pos)
        #println(i)
        cell = equidistant_normal_points_matrix(cell, 1.0, 18.0)
        #CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=3, color=KubaPhD.ColorsCBfriendly[:Black], joinstyle=:round)
        CairoMakie.poly!(ax, cell[1,:], cell[2,:], color=:black, strokecolor=:black, strokewidth=2, joinstyle=:round, stroke_depth_shift = +1e-3)
        CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=3, color=:black, joinstyle=:round)
    end
end

function PlotWallEmbedDensities(sol, embedded_pos, ξ_all, Exp_Ot)

    fig1 = Figure(size=(600, 600))
    ga = fig1[1,1] = GridLayout()

    WallDensityExample_ax = CairoMakie.Axis(ga[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", xlabelsize = 22, ylabelsize = 22, aspect = 1, xticklabelsize = 18, yticklabelsize = 18, yticks = (-300:100:300, ["-300", "-200", "-100", "0", "100", "200", "300"]))
    #colors = [:dodgerblue, :darkorange2, :green, :mediumorchid2]
    PlotOsteonSim!(WallDensityExample_ax, sol, embedded_pos)
    hidedecorations!(WallDensityExample_ax;grid=false)


    fig2 = Figure(size=(600, 600))
    colors = reverse(greyscale_colors(size(sol.t, 1)))
    ga = fig2[1,1] = GridLayout()
    data_ax = CairoMakie.Axis(ga[1,1], xlabelsize = 22, ylabelsize = 22, xticklabelsize = 18, yticklabelsize = 18, aspect = 1, limits=(0.5,7.5,-0.0001, 0.0021), xticks = (1.0:1.0:7.0, ["T", "2T", "3T", "4T", "5T", "6T", "7T"]))
    for ii in 1:size(ξ_all[:,1],1)
        CairoMakie.boxplot!(data_ax, ii*ones(size(ξ_all[1,:],1)),ξ_all[ii,:], color=colors[ii], whiskerwidth=0.5, show_outliers=false, strokecolor = :black, strokewidth = 1.5)
    end
    CairoMakie.lines!(data_ax, 0:size(ξ_all[:,1],1)+1, ones(size(0:size(ξ_all[:,1],1)+1)).*Exp_Ot, color=:red, linestyle=:dash, linewidth=3, label=L"\text{Target Ot}")
    for ii in 1:size(ξ_all[:,1],1)
        if ii == 1
            CairoMakie.scatter!(data_ax, 1, mean(ξ_all[1,:]), color=:dodgerblue, marker=:star5, strokecolor = :black, strokewidth = 1, markersize=20, label=L"\text{Average density}")
        else
            CairoMakie.scatter!(data_ax, ii, mean(ξ_all[ii,:]), color=:dodgerblue, marker=:star5, strokecolor = :black, strokewidth = 1, markersize=20)
        end
    end

    CairoMakie.axislegend(data_ax, position = :rt,  labelsize = 24)
    #CairoMakie.scatter!(data_ax, 0:size(ξ_all[:,1],1)+1,ones(size(0:size(ξ_all[:,1],1)+1)).*mean(total_ξ), color=KubaPhD.BlueRedColors[:MidBlue2], marker=:diamond, markersize=20, label="Total density")

    return fig1, fig2
end


function PlotWThAsymmetry(WTh_results, all_solutions, all_embedded_cell_pos)

    smodel = "box"

    if smodel == "linear"
        formula_small = StatsModels.@formula(ξ_small ~ 1 + WThRatio)
        formula_large = StatsModels.@formula(ξ_large ~ 1 + WThRatio)
        model_small = GLM.lm(formula_small, WTh_results)
        model_large = GLM.lm(formula_large, WTh_results)

        # Plotting

        fig1 = Figure(size=(700, 700))
        ga = fig1[1, 1] = GridLayout()

        ax = CairoMakie.Axis(
            ga[1, 1],
            xlabel=L"$\text{W.Th asymmetry ratio}$",
            ylabel=L"$\xi$ [cells/μm$^2$]",
            xlabelsize=22,
            ylabelsize=22,
            xticklabelsize=18,
            yticklabelsize=18,
            aspect=1,
            limits = (0.25, 0.95, -0.0003, 0.004)
        )

        CairoMakie.scatter!(ax, WTh_results.WThRatio, WTh_results.ξ_small, markersize=10, color=:blue, alpha=0.5, marker=:circle, label="Small Sector")
        CairoMakie.scatter!(ax, WTh_results.WThRatio, WTh_results.ξ_large, markersize=10, color=:red, alpha=0.5, marker=:circle, label="Large Sector")
        CairoMakie.lines!(ax, WTh_results.WThRatio, vcat(predict(model_small, WTh_results)...), color=:blue, linewidth=4, label="LS-Fit Small")
        CairoMakie.lines!(ax, WTh_results.WThRatio, vcat(predict(model_large, WTh_results)...), color=:red, linewidth=4, label="LS-Fit Large")
        CairoMakie.axislegend(ax, position=:rt, fontsize=12)
    else
        fig1 = Figure(size=(900, 900))
        ga = fig1[1, 1] = GridLayout()
        ax = CairoMakie.Axis(
            ga[1, 1],
            #xlabel=L"$\text{W.Th asymmetry ratio}$",
            #ylabel=L"$\xi$ [cells/μm$^2$]",
            xlabelsize=22,
            ylabelsize=22,
            xticklabelsize=24,
            yticklabelsize=24,
            aspect=1,
            limits = (0.3, 1.0, -0.0003, 0.0045),
            xticks = ([0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95], ["[0.3,0.4)", "[0.4,0.5)", "[0.5,0.6)", "[0.6,0.7)","[0.7,0.8)", "[0.8,0.9)", "[0.9,1.0)"]),
        )

        # Define the ranges for grouping
        ranges = 0.2:0.1:1.0

        # Create a new column in the WTh_results_table for the group each WThRatio belongs to
        function custom_cut(values, ranges)
            groups = Vector{Int}(undef, length(values))
            for i in 1:length(values)
                for j in 1:(length(ranges) - 1)
                    if ranges[j] <= values[i] < ranges[j + 1]
                        groups[i] = j
                        break
                    end
                end
            end
            return groups
        end

        WTh_results[!, :WThRatio_group] = custom_cut(WTh_results.WThRatio, collect(ranges))

        # Create a separate DataFrame for each group
        grouped_dataframes = Dict{Int, DataFrame}()
        for group in unique(WTh_results.WThRatio_group)
            grouped_dataframes[group] = filter(row -> row.WThRatio_group == group, WTh_results)
        end
        # Plot the grouped data
        bins = [0.225, 0.325, 0.425, 0.525, 0.625, 0.725, 0.825, 0.925, 0.275, 0.375, 0.475, 0.575, 0.675, 0.775, 0.875, 0.975]
        #CairoMakie.boxplot!(ax, bins[1]*ones(size(grouped_dataframes[1],1)), grouped_dataframes[1].ξ_small, color=:blue, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[2]*ones(size(grouped_dataframes[2],1)), grouped_dataframes[2].ξ_small, color=:blue, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[3]*ones(size(grouped_dataframes[3],1)), grouped_dataframes[3].ξ_small, color=:blue, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[4]*ones(size(grouped_dataframes[4],1)), grouped_dataframes[4].ξ_small, color=:blue, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[5]*ones(size(grouped_dataframes[5],1)), grouped_dataframes[5].ξ_small, color=:blue, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[6]*ones(size(grouped_dataframes[6],1)), grouped_dataframes[6].ξ_small, color=:blue, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[7]*ones(size(grouped_dataframes[7],1)), grouped_dataframes[7].ξ_small, color=:blue, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[8]*ones(size(grouped_dataframes[8],1)), grouped_dataframes[8].ξ_small, color=:blue, whiskerwidth=0.5, width=0.05, show_outliers=false)


        #CairoMakie.boxplot!(ax, bins[9]*ones(size(grouped_dataframes[1],1)), grouped_dataframes[1].ξ_large, color=:red, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[10]*ones(size(grouped_dataframes[2],1)), grouped_dataframes[2].ξ_large, color=:red, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[11]*ones(size(grouped_dataframes[3],1)), grouped_dataframes[3].ξ_large, color=:red, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[12]*ones(size(grouped_dataframes[4],1)), grouped_dataframes[4].ξ_large, color=:red, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[13]*ones(size(grouped_dataframes[5],1)), grouped_dataframes[5].ξ_large, color=:red, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[14]*ones(size(grouped_dataframes[6],1)), grouped_dataframes[6].ξ_large, color=:red, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[15]*ones(size(grouped_dataframes[7],1)), grouped_dataframes[7].ξ_large, color=:red, whiskerwidth=0.5, width=0.05, show_outliers=false)
        CairoMakie.boxplot!(ax, bins[16]*ones(size(grouped_dataframes[8],1)), grouped_dataframes[8].ξ_large, color=:red, whiskerwidth=0.5, width=0.05, show_outliers=false)

    end
    hidedecorations!(ax; grid=false)


    # Plotting growth simulation results
    fig2 = Figure(size=(1400, 700))
    ga = fig2[1,1] = GridLayout()
    gb = fig2[1,2] = GridLayout()
    gc = fig2[1,3] = GridLayout()
    gd = fig2[1,4] = GridLayout()
    ge = fig2[2,1] = GridLayout()
    gf = fig2[2,2] = GridLayout()
    gg = fig2[2,3] = GridLayout()
    gh = fig2[2,4] = GridLayout()

    growth_sim_ax_1 = CairoMakie.Axis(ga[1,1], ylabel=L"$y$ [μm]", xlabelsize = 22, ylabelsize = 22, aspect = 1,xticklabelsize=18, yticklabelsize=18)
    growth_sim_ax_2 = CairoMakie.Axis(gb[1,1], xlabelsize = 22, ylabelsize = 22, aspect = 1,xticklabelsize=18, yticklabelsize=18)
    growth_sim_ax_3 = CairoMakie.Axis(gc[1,1], xlabelsize = 22, ylabelsize = 22, aspect = 1,xticklabelsize=18, yticklabelsize=18)
    growth_sim_ax_4 = CairoMakie.Axis(gd[1,1], xlabelsize = 22, ylabelsize = 22, aspect = 1,xticklabelsize=18, yticklabelsize=18)
    growth_sim_ax_5 = CairoMakie.Axis(ge[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", xlabelsize = 22, ylabelsize = 22, aspect = 1,xticklabelsize=18, yticklabelsize=18)
    growth_sim_ax_6 = CairoMakie.Axis(gf[1,1], xlabel=L"$x$ [μm]", xlabelsize = 22, ylabelsize = 22, aspect = 1,xticklabelsize=18, yticklabelsize=18)
    growth_sim_ax_7 = CairoMakie.Axis(gg[1,1], xlabel=L"$x$ [μm]", xlabelsize = 22, ylabelsize = 22, aspect = 1,xticklabelsize=18, yticklabelsize=18)
    growth_sim_ax_8 = CairoMakie.Axis(gh[1,1], xlabel=L"$x$ [μm]", xlabelsize = 22, ylabelsize = 22, aspect = 1,xticklabelsize=18, yticklabelsize=18)

    all_axes = [growth_sim_ax_1, growth_sim_ax_2, growth_sim_ax_3, growth_sim_ax_4, growth_sim_ax_5, growth_sim_ax_6, growth_sim_ax_7, growth_sim_ax_8]
    idxes = [250, 45, 1, 15, 53, 172, 60, 30]#[249, 24, 1, 194, 230, 130, 32, 255]#[199, 225, 210, 157, 186, 59, 243, 240] #[33, 28, 25, 15, 11, 31, 24, 50] # selected simulations with WTh ratios

    all_sols = all_solutions;
    all_embedded_pos = all_embedded_cell_pos;

    for ii in eachindex(all_axes)
        PlotOsteonSim!(all_axes[ii], all_sols[idxes[ii]], all_embedded_pos[idxes[ii]])
        hidedecorations!(all_axes[ii]; grid=false)
    end

    return fig1, fig2
end