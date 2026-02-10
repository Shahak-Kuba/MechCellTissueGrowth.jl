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
end