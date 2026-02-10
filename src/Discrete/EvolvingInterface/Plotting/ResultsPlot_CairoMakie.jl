

function PlotGrowthSimulation(data, embedded_cells)
    # Figure 1: Growth simulation, number of cells, embedded cells, interface density, filled tissue area
    fig1 = Figure(size=(700, 700))
    ga = fig1[1,1] = GridLayout()
    growth_sim_ax = Axis(ga[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", title="Growth simulation", xlabelsize = 18, ylabelsize = 18, aspect=1)
    plotGrowthSimulation_no_cb!(growth_sim_ax, data.u, data.Density, embedded_cells, :black, :black)

    return fig1
end

function PlotGrowthSimulation_ImageStack(data, embedded_cells)
    fig1 = Figure(size=(700, 700), backgroundcolor = RGBA(0, 0, 0, 0.1))
    fig2 = Figure(size=(700, 700), backgroundcolor = RGBA(0, 0, 0, 0.1))
    
    ga = fig1[1,1] = GridLayout()
    gb = fig2[1,1] = GridLayout()
    axisLimits = maximum([maximum(data.u[1][1,:]), maximum(data.u[1][2,:])]) + 20.0

    Interface_ax = Axis(ga[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", xlabelsize = 18, ylabelsize = 18, aspect=1, limits=(-axisLimits,axisLimits,-axisLimits,axisLimits), backgroundcolor = RGBA(0, 0, 0, 0))
    Embedded_cells_ax = Axis(gb[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", xlabelsize = 18, ylabelsize = 18, aspect=1, limits=(-axisLimits,axisLimits,-axisLimits,axisLimits), backgroundcolor = RGBA(0, 0, 0, 0))

    plotInterface_no_cb!(Interface_ax, data.u, data.Density,1, :blue)
    plotInterface_no_cb!(Interface_ax, data.u, data.Density,size(data.u,1), :blue)
    
    
    centroid = KubaPhD.calc_shape_centroid(data, size(data.u,1))
    scatter!(Interface_ax, centroid[1], centroid[2], markersize=15, color=:green)

    for cell in embedded_cells
        cell = equidistant_normal_points_matrix(cell, 1.0)
        CairoMakie.lines!(Embedded_cells_ax, cell[1,:], cell[2,:], linewidth=3, color=:red)
        CairoMakie.lines!(Embedded_cells_ax, cell[1,:], cell[2,:], linewidth=2, color=:red)
    end

    return fig1, fig2
end

function PlotSimulationOutputs_Cell_Culture(data, embedded_cells, embedded_cell_count, embed_rates, Embed, CellMech)

    # Figure 1: Growth simulation, number of cells, embedded cells, interface density, filled tissue area
    fig1 = Figure(size=(700, 600))
    ga = fig1[1,1] = GridLayout()
    gb = fig1[1,2] = GridLayout()
    gc = fig1[2,1] = GridLayout()
    gd = fig1[2,2] = GridLayout()

    growth_sim_ax = Axis(ga[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", title="Growth simulation", xlabelsize = 18, ylabelsize = 18)
    cell_count_ax = Axis(gb[1,1], xlabel=L"$t$ [days]", ylabel=L"N", title="Number of cells", yticklabelcolor = :blue, xlabelsize = 18, ylabelsize = 18)
    embedded_cell_count_ax = Axis(gb[1,1], yaxisposition=:right, ylabel=L"N_{embed}", yticklabelcolor = :red, xlabelsize = 18, ylabelsize = 18)
    interface_density_ax = Axis(gc[1,1], xlabel=L"$t$ [days]", ylabel=L"$q$ [1/μm]", title="Mean interface density", xlabelsize = 18, ylabelsize = 18)
    tissue_area_ax = Axis(gd[1,1], xlabel=L"$t$ [days]", ylabel=L"$\Omega$ [μm²]", title="Filled tissue area over time", xlabelsize = 18, ylabelsize = 18)

    plotGrowthSimulation_no_cb!(growth_sim_ax, data.u, data.Density, embedded_cells)

    plotCellCount!(cell_count_ax, embedded_cell_count_ax, data.t, data.CellCount, embedded_cell_count[1])

    plotInterfaceDensity!(interface_density_ax, data.t, data.Density)

    plotTissueArea!(tissue_area_ax, data.t, data.Ω)
    #plotEmbedmentRate!(tissue_area_ax, data.t, Embed.rate, embed_rates[1], CellMech.kf)

    # adding axes labels
    for (label, layout) in zip(["a)", "b)", "c)", "d)"], [ga, gb, gc, gd])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 18,
            padding = (0, 30, 15, 0),
            halign = :right)
    end


    # Figure 2: Cell orientation plot

    fig2 = Figure(size=(900, 250))
    ga2 = fig2[1,1] = GridLayout()
    gb2 = fig2[1,2] = GridLayout()
    gc2 = fig2[1,3] = GridLayout()

    idx1 = 10; t1 = floor(data.t[idx1])
    idx2 = Int(round(length(data.t)/2)); t2 = floor(data.t[idx2])
    idx3 = length(data.t); t3 = floor(data.t[end])

    cell_orientation_ax_1 = PolarAxis(ga2[1,1], thetalimits = (0, pi), title =L"$t = %$t1$ days", thetaticks = (0:π/4:π, ["0°", "45°", "90°", "135°", "180°"]))
    cell_orientation_ax_2 = PolarAxis(gb2[1,1], thetalimits = (0, pi), title =L"$t = %$t2$ days", thetaticks = (0:π/4:π, ["0°", "45°", "90°", "135°", "180°"]))
    cell_orientation_ax_3 = PolarAxis(gc2[1,1], thetalimits = (0, pi), title =L"$t = %$t3$ days", thetaticks = (0:π/4:π, ["0°", "45°", "90°", "135°", "180°"]))
    supertitle = Label(fig2[0, :], L"Orientation of embedded cells at different values of $t$", fontsize = 24)

    plotCellOrientation_at_t!(cell_orientation_ax_1, embedded_cells, embedded_cell_count, idx1)
    plotCellOrientation_at_t!(cell_orientation_ax_2, embedded_cells, embedded_cell_count, idx2)
    plotCellOrientation_at_t!(cell_orientation_ax_3, embedded_cells, embedded_cell_count, idx3)

    return fig1, fig2 

end


function PlotSimulationOutputs_Bone(data, embedded_cells, embedded_cell_count, embed_rates, CellMech, Embed)

    # Figure 1: Growth simulation, number of cells, embedded cells, interface density, filled tissue area
    fig1 = Figure(size=(700, 600))
    ga = fig1[1,1] = GridLayout()
    gb = fig1[1,2] = GridLayout()
    gc = fig1[2,1] = GridLayout()
    gd = fig1[2,2] = GridLayout()

    growth_sim_ax = Axis(ga[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", title="Growth simulation", xlabelsize = 18, ylabelsize = 18)
    cell_count_ax = Axis(gc[1,1], xlabel=L"$t$ [days]", ylabel=L"N", title="Number of cells", yticklabelcolor = :blue, xlabelsize = 18, ylabelsize = 18)
    embedded_cell_count_ax = Axis(gc[1,1], yaxisposition=:right, ylabel=L"N_{embed}", yticklabelcolor = :red, xlabelsize = 18, ylabelsize = 18)
    Ot_vs_simOt_ax = Axis(gb[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", title="Embed distribution heatmap", xlabelsize = 18, ylabelsize = 18)
    Ot_rate_ax = Axis(gd[1,1], xlabel=L"$t$ [days]", ylabel=L"$⟨E\;⟩/k$_{f}", title="Mean embedment rate", xlabelsize = 18, ylabelsize = 18)

    plotGrowthSimulation_no_cb!(growth_sim_ax, data.u, data.Density, embedded_cells)

    plotEmbedDensity!(Ot_vs_simOt_ax, data.t, embedded_cell_count, data.Ω, CellMech, Embed)

    plotCellCount!(cell_count_ax, embedded_cell_count_ax, data.t, data.CellCount, embedded_cell_count[1])

    #plotTissueArea!(tissue_area_ax, data.t, data.Ω)
    plotEmbedmentRate!(Ot_rate_ax, data.t, Embed.rate, embed_rates[1], CellMech.kf)

    # adding axes labels
    for (label, layout) in zip(["a)", "b)", "c)", "d)"], [ga, gb, gc, gd])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 18,
            padding = (0, 30, 15, 0),
            halign = :right)
    end

    return fig1

end

function PlotDiscreteVsContinuum(all_disc_outputs, all_disc_embed_pos, continuum_output, cmap, clims)

    fig1 = Figure(size=(900, 800))
    ga = fig1[1,1] = GridLayout()
    gb = fig1[1,2] = GridLayout()
    gc = fig1[1,3] = GridLayout()
    gd = fig1[2,1] = GridLayout()
    ge = fig1[2,2] = GridLayout() # Continuum model plot
    gf = fig1[2,3] = GridLayout()
    gg = fig1[3,:] = GridLayout()
    
    disc_1_ax = Axis(ga[1,1], ylabel=L"$y$ [μm]", title="Realisation 1", xlabelsize = 18, ylabelsize = 18)
    disc_2_ax = Axis(gb[1,1], title="Realisation 2", xlabelsize = 18, ylabelsize = 18)
    disc_3_ax = Axis(gc[1,1], title="Realisation 3", xlabelsize = 18, ylabelsize = 18)
    disc_4_ax = Axis(gd[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", title="Realisation 4", xlabelsize = 18, ylabelsize = 18)
    cont_ax = Axis(gf[1,1], xlabel=L"$x$ [μm]", title="Continuum simulation", xlabelsize = 18, ylabelsize = 18)
    disc_5_ax = Axis(ge[1,1], xlabel=L"$x$ [μm]", title="Realisation 5", xlabelsize = 18, ylabelsize = 18)
    mean_int_density_ax = Axis(gg[1,1], xlabel=L"$t$ [days]", ylabel=L"$\bar{q}$ [$1/\mu$m]", title="Mean interface density", xlabelsize = 18, ylabelsize = 18)


    plotGrowthSimulation!(disc_1_ax, all_disc_outputs[1].u, all_disc_outputs[1].Density, all_disc_embed_pos[1], (200.0, 200.0), clims, cmap)
    plotGrowthSimulation!(disc_2_ax, all_disc_outputs[2].u, all_disc_outputs[2].Density, all_disc_embed_pos[2], (200.0, 200.0), clims, cmap)
    plotGrowthSimulation!(disc_3_ax, all_disc_outputs[3].u, all_disc_outputs[3].Density, all_disc_embed_pos[3], (200.0, 200.0), clims, cmap)
    plotGrowthSimulation!(disc_4_ax, all_disc_outputs[4].u, all_disc_outputs[4].Density, all_disc_embed_pos[4], (200.0, 200.0), clims, cmap) 
    plotGrowthSimulation!(disc_5_ax, all_disc_outputs[5].u, all_disc_outputs[5].Density, all_disc_embed_pos[5], (200.0, 200.0), clims, cmap)
    plotContinuumSimulation!(cont_ax, continuum_output.θ, continuum_output.R, continuum_output.ρ, cmap, clims)
    CairoMakie.Colorbar(fig1[1:2,4], limits = clims, colormap = cmap, label = L"Cell density $q$ [$1/\mu$m]", labelsize = 20 ,labelpadding = 10, ticklabelsize = 12, width = 20)

    plotMeanInterfaceDensity_DiscVsCont!(mean_int_density_ax, all_disc_outputs, continuum_output)
    CairoMakie.Legend(fig1[3,4], mean_int_density_ax)
    

    # adding axes labels
    for (label, layout) in zip(["a)", "b)", "c)", "d)", "e)", "f)", "g)"], [ga, gb, gc, gd, ge, gf, gg])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 18,
            padding = (0, 30, 15, 0),
            halign = :right)
    end

    return fig1
end

function PlotDiscreteVsContinuum(all_disc_outputs_m_1, all_disc_embed_pos_1, continuum_output, cmap, clims, all_disc_outputs_m_2, all_disc_outputs_m_3)

    # Figure 1: Growth simulation, number of cells, embedded cells, interface density, filled tissue area
    fig1 = Figure(size=(900, 1000))
    ga = fig1[1,1] = GridLayout()
    gb = fig1[1,2] = GridLayout()
    gc = fig1[1,3] = GridLayout()
    gd = fig1[2,1] = GridLayout()
    ge = fig1[2,2] = GridLayout() # Continuum model plot
    gf = fig1[2,3] = GridLayout()
    gg = fig1[3,:] = GridLayout()
    gh = fig1[4,:] = GridLayout()

    disc_1_ax = Axis(ga[1,1], ylabel=L"$y$ [μm]", title="Realisation 1", xlabelsize = 18, ylabelsize = 18)
    disc_2_ax = Axis(gb[1,1], title="Realisation 2", xlabelsize = 18, ylabelsize = 18)
    disc_3_ax = Axis(gc[1,1], title="Realisation 3", xlabelsize = 18, ylabelsize = 18)
    disc_4_ax = Axis(gd[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", title="Realisation 4", xlabelsize = 18, ylabelsize = 18)
    cont_ax = Axis(gf[1,1], xlabel=L"$x$ [μm]", title="Continuum simulation", xlabelsize = 18, ylabelsize = 18)
    disc_5_ax = Axis(ge[1,1], xlabel=L"$x$ [μm]", title="Realisation 5", xlabelsize = 18, ylabelsize = 18)
    mean_int_density_ax = Axis(gg[1,1], xlabel=L"$t$ [days]", ylabel=L"$\bar{q}$ [$1/\mu$m]", title="Mean interface density", xlabelsize = 18, ylabelsize = 18)
    mean_area_ax = Axis(gh[1,1], xlabel=L"$t$ [days]", ylabel=L"$\Omega$ [$\mu$m²]", title="Area of tissue generated", xlabelsize = 18, ylabelsize = 18)


    plotGrowthSimulation_3_timepoints!(disc_1_ax, all_disc_outputs_m_1[1].u, all_disc_outputs_m_1[1].Density, all_disc_embed_pos_1[1], (200.0, 200.0), clims, cmap)
    plotGrowthSimulation_3_timepoints!(disc_2_ax, all_disc_outputs_m_1[2].u, all_disc_outputs_m_1[2].Density, all_disc_embed_pos_1[2], (200.0, 200.0), clims, cmap)
    plotGrowthSimulation_3_timepoints!(disc_3_ax, all_disc_outputs_m_1[3].u, all_disc_outputs_m_1[3].Density, all_disc_embed_pos_1[3], (200.0, 200.0), clims, cmap)
    plotGrowthSimulation_3_timepoints!(disc_4_ax, all_disc_outputs_m_1[4].u, all_disc_outputs_m_1[4].Density, all_disc_embed_pos_1[4], (200.0, 200.0), clims, cmap) 
    plotGrowthSimulation_3_timepoints!(disc_5_ax, all_disc_outputs_m_1[5].u, all_disc_outputs_m_1[5].Density, all_disc_embed_pos_1[5], (200.0, 200.0), clims, cmap)
    plotContinuumSimulation!(cont_ax, continuum_output.θ, continuum_output.R, continuum_output.ρ, cmap, clims)
    CairoMakie.Colorbar(fig1[1:2,4], limits = clims, colormap = cmap, label = L"Cell density $q$ [$1/\mu$m]", labelsize = 20 ,labelpadding = 10, ticklabelsize = 12, width = 20)

    #plotMeanInterfaceDensity_DiscVsCont!(mean_int_density_ax, all_disc_outputs, continuum_output)
    plotMeanInterfaceDensity_DiscVsCont!(mean_int_density_ax, all_disc_outputs_m_2, all_disc_outputs_m_3, all_disc_outputs_m_1, continuum_output)
    CairoMakie.Legend(fig1[3:4,4], mean_int_density_ax)

    #plotTissueArea_DiscVsCont!(mean_area_ax, all_disc_outputs_m_1, continuum_output)
    plotTissueArea_DiscVsCont!(mean_area_ax, all_disc_outputs_m_2, all_disc_outputs_m_3, all_disc_outputs_m_1, continuum_output)
    

    # adding axes labels
    for (label, layout) in zip(["a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"], [ga, gb, gc, gd, ge, gf, gg, gh])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 18,
            padding = (0, 30, 15, 0),
            halign = :right)
    end

    return fig1
end