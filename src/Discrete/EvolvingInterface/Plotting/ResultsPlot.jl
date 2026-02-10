

function PlotSimulationOutputs_Cell_Culture(data, embedded_cells, embedded_cell_count, embed_rates)

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

    plotGrowthSimulation_no_cb!(growth_sim_ax, data.u, data.Density, embedded_cells, (250.0, 250.0), (0.05,0.12), :dense)

    plotCellCount!(cell_count_ax, embedded_cell_count_ax, data.t, data.CellCount, embedded_cell_count[1])

    plotInterfaceDensity!(interface_density_ax, data.t, data.Density)

    #plotTissueArea!(tissue_area_ax, data.t, data.Ω)
    plotEmbedmentRate!(tissue_area_ax, data.t, 0.5, embed_rates[1], 87.842)

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
    Ot_heat_map_ax = Axis(gb[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", title="Embed distribution heatmap", xlabelsize = 18, ylabelsize = 18)
    Ot_rate_ax = Axis(gd[1,1], xlabel=L"$t$ [days]", ylabel=L"$⟨E\;⟩/k$_{f}", title="Mean embedment rate", xlabelsize = 18, ylabelsize = 18)

    plotGrowthSimulation_no_cb!(growth_sim_ax, data.u, data.Density, embedded_cells, (250.0, 250.0), (0.05,0.12), :dense)

    tr = plotEmbedCellDensityMap!(Ot_heat_map_ax, data)
    Colorbar(fig1[1, 2][1, 2], tr, label = "Embedded cell count", labelpadding = 10, ticklabelsize = 12)

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