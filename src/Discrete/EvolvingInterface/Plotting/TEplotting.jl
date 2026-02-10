
function blackscale_colors(N::Integer)
    N ≥ 1 || throw(ArgumentError("N must be ≥ 1"))
    vals = range(0.0, 0.5; length = N)
    return [RGB(v, v, v) for v in vals]
end

function PlotTE_Simulation(sol, embedded_pos)
    fig = Figure(size=(800,800))
    ga = fig[1,1] = GridLayout()
    ax = CairoMakie.Axis(ga[1,1], xlabel=L"$x$ [μm]", ylabel=L"$y$ [μm]", xlabelsize = 22, xlabelvisible = false, ylabelvisible= false, ylabelsize = 22, aspect = 1, xticklabelsize = 18, yticklabelsize = 18, xticklabelsvisible=false, yticklabelsvisible=false, backgroundcolor=:white)
    
    TissueColor = KubaPhD.ColorsCBfriendly[:Olive]
    InterfaceColor1 = RGB(38/255, 90/255, 166/255)
    InterfaceColor2 = RGB(90/255, 150/255, 90/255)
    InterfaceColor3 = RGB(245/255, 158/255, 11/255)
    CellColor = RGB(220/255, 50/255, 47/255)
    CentreFillColor = RGB(1, 1, 1)    


    fake_scaffold = false

    if fake_scaffold
        # find min and max x and y values to make a fake scaffold
        min_x = minimum(sol.u[1][1,:])
        max_x = maximum(sol.u[1][1,:])
        min_y = minimum(sol.u[1][2,:])
        max_y = maximum(sol.u[1][2,:])
        scaffold_width = 20.0
        # make a fake scaffold around the initial tissue
        scaffold_left_x = [min_x - scaffold_width, min_x, min_x, min_x - scaffold_width]
        scaffold_left_y = [min_y - 3*scaffold_width, min_y - 3*scaffold_width, max_y + 3*scaffold_width, max_y + 3*scaffold_width]
        CairoMakie.poly!(ax, scaffold_left_x, scaffold_left_y, color=KubaPhD.ColorsCBfriendly[:NeutralGray], strokecolor=:black, strokewidth=2, joinstyle=:round, stroke_depth_shift = +1e-3)
        scaffold_right_x = [max_x, max_x + scaffold_width, max_x + scaffold_width, max_x]
        scaffold_right_y = [min_y - 3*scaffold_width, min_y - 3*scaffold_width, max_y + 3*scaffold_width, max_y + 3*scaffold_width]
        CairoMakie.poly!(ax, scaffold_right_x, scaffold_right_y, color=KubaPhD.ColorsCBfriendly[:NeutralGray], strokecolor=:black, strokewidth=2, joinstyle=:round, stroke_depth_shift = +1e-3)
        scaffold_top_x = [min_x - 3*scaffold_width, max_x + 3*scaffold_width, max_x + 3*scaffold_width, min_x - 3*scaffold_width]
        scaffold_top_y = [max_y, max_y, max_y + scaffold_width, max_y + scaffold_width]
        CairoMakie.poly!(ax, scaffold_top_x, scaffold_top_y, color=KubaPhD.ColorsCBfriendly[:NeutralGray], strokecolor=:black, strokewidth=2, joinstyle=:round, stroke_depth_shift = +1e-3)
        scaffold_bottom_x = [min_x - 3*scaffold_width, max_x + 3*scaffold_width, max_x + 3*scaffold_width, min_x - 3*scaffold_width]
        scaffold_bottom_y = [min_y - scaffold_width, min_y - scaffold_width, min_y, min_y]
        CairoMakie.poly!(ax, scaffold_bottom_x, scaffold_bottom_y, color=KubaPhD.ColorsCBfriendly[:NeutralGray], strokecolor=:black, strokewidth=2, joinstyle=:round, stroke_depth_shift = +1e-3)
        limits!(ax, min_x - 3*scaffold_width, max_x + 3*scaffold_width, min_y - 3*scaffold_width, max_y + 3*scaffold_width)
    end


    poly!(ax, KubaPhD.convert_coordinates_to_tuples(sol.u[1]), color=CentreFillColor, strokecolor = :black, strokewidth = 2, joinstyle=:round, alpha=0.3, stroke_depth_shift = +1e-3)
    #poly!(ax, KubaPhD.convert_coordinates_to_tuples(sol.u[4]), color=TissueColor, strokecolor = :black, strokewidth = 2, joinstyle=:round, alpha=0.3, stroke_depth_shift = +1e-3)
    #poly!(ax, KubaPhD.convert_coordinates_to_tuples(sol.u[8]), color=TissueColor, strokecolor = :black, strokewidth = 2, joinstyle=:round, alpha=0.3, stroke_depth_shift = +1e-3)
    poly!(ax, KubaPhD.convert_coordinates_to_tuples(sol.u[end]), color=CentreFillColor, strokecolor = :black, strokewidth = 2, joinstyle=:round,  stroke_depth_shift = +1e-3)
    #for ii in axes(sol.t, 1)
    #   KubaPhD.plotInterface_no_cb!(ax,sol.u,sol.Density,ii,InterfaceColor)
    #end

    for cell in embedded_pos
        cell = KubaPhD.equidistant_normal_points_matrix(cell, 1.5, 15.0)
        #CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=3, color=KubaPhD.ColorsCBfriendly[:Black], joinstyle=:round)
        CairoMakie.poly!(ax, cell[1,:], cell[2,:], color=CellColor, strokecolor=:black, strokewidth=1, joinstyle=:round, stroke_depth_shift = +1e-3, alpha=0.8)
        CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=1, color=:black, joinstyle=:round)
    end

    KubaPhD.plotInterface_no_cb!(ax,sol.u,sol.Density,1,:black)
    KubaPhD.plotInterface_no_cb!(ax,sol.u,sol.Density,5,InterfaceColor1)
    KubaPhD.plotInterface_no_cb!(ax,sol.u,sol.Density,8,InterfaceColor2)
    KubaPhD.plotInterface_no_cb!(ax,sol.u,sol.Density,length(sol.u),InterfaceColor3)
    return fig
end