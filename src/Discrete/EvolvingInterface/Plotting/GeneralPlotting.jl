

const myColors = Dict(
    :DarkBlue => RGB(0.23921568627450981, 0.35294117647058826, 0.5019607843137255),
    :DeepNavy => RGB(0.06274509803921569, 0.08627450980392157, 0.30980392156862746),
    :SlateBlue => RGB(0.41568627450980394, 0.35294117647058826, 0.803921568627451),
    :SkyBlue => RGB(0.3411764705882353, 0.6078431372549019, 0.7568627450980392),
    :MidBlue => RGB(0.49411764705882355, 0.6588235294117647, 0.7450980392156863),
    :LightBlue => RGB(0.596078431372549, 0.7568627450980392, 0.8509803921568627),
    :LightBlue1 => RGB(0.5137254901960784, 0.7098039215686275, 0.8196078431372549),
    :DarkRed => RGB(0.42745098039215684, 0.12549019607843137, 0.1450980392156863),
    :Red => RGB(0.5019607843137255, 0.0, 0.0),
    :Orange => RGB(0.9058823529411765, 0.5647058823529412, 0.29411764705882354),
    :BurntOrange => RGB(0.8509803921568627, 0.37254901960784315, 0.00784313725490196),
    :DarkGrey => RGB(0.1607843137254902, 0.19607843137254902, 0.2549019607843137),
    :PeachYellow => RGB(0.9294117647058824, 0.8117647058823529, 0.5568627450980392),
    :Turquoise => RGB(0.40784313725490196, 0.6901960784313725, 0.6705882352941176),
    :YellowSaffron => RGB(0.9725490196078431, 0.7764705882352941, 0.18823529411764706),
    :Magenta => RGB(0.9058823529411765, 0.1607843137254902, 0.5411764705882353)
)

const ColorsCBfriendly = Dict(
    :DodgerBlue => RGB(0.11764705882352941, 0.56470588235294117, 1.00000000000000000),  # #1E90FF
    :DeepBlue => RGB(0.00000000000000000, 0.44705882352941179, 0.69803921568627447),  # #0072B2
    :Orange => RGB(0.90196078431372551, 0.62352941176470589, 0.00000000000000000),  # #E69F00
    :SkyBlue => RGB(0.33725490196078434, 0.70588235294117652, 0.91372549019607840),  # #56B4E9
    :BluishGreen => RGB(0.00000000000000000, 0.61960784313725492, 0.45098039215686275),  # #009E73
    :Olive => RGB(0.47843137254901963, 0.60392156862745094, 0.00392156862745098),  # #7A9A01
    :Vermillion => RGB(0.83529411764705885, 0.36862745098039218, 0.00000000000000000),  # #D55E00
    :ReddishPurple => RGB(0.80000000000000004, 0.47450980392156861, 0.65490196078431373),  # #CC79A7
    :NeutralGray => RGB(0.49803921568627452, 0.49803921568627452, 0.49803921568627452),  # #7F7F7F
    :Black => RGB(0.00000000000000000, 0.00000000000000000, 0.00000000000000000)  # #000000
)



function data_for_plotting(u,var,index)
    xp = [u[index][1, :]; u[index][1,1]]
    yp = [u[index][2, :]; u[index][2,1]]
    if typeof(var) == Vector{Vector{Float64}}
        varp = [var[index]; var[index][1]]
    else
        varp = [var[index]; var[index][1]].data
    end
    return xp, yp, varp
end

function plotInterface!(ax, u,var,index,cbar_limits,cmap)
    x,y,z = data_for_plotting(u, var, index)
    CairoMakie.lines!(ax,x,y,color=z,linewidth=4, colorrange=cbar_limits, colormap=cmap)
end

function plotGrowthSimulation!(ax, u, var, embedded_cell_pos, axis_limits, cbar_limits, cmap)

    # plotting the initial interface
    plotInterface!(ax, u,var,1,cbar_limits,cmap)

    # plotting the embedded cells
    for cell in embedded_cell_pos
        cell = equidistant_normal_points_matrix(cell, 1.0)
        CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=3, color=:black)
        CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=2, color=:red)
    end

    # plotting the final interface
    plotInterface!(ax, u,var,length(u),cbar_limits,cmap)


end


function plotInterface_no_cb!(ax, u,var,index,linecolor)
    x,y,z = data_for_plotting(u, var, index)
    CairoMakie.lines!(ax,x,y,color=linecolor,linewidth=5)
end


function plotGrowthSimulation_no_cb!(ax, u, var, embedded_cell_pos, clr1, clr2)

    # plotting the initial interface
    plotInterface_no_cb!(ax,u,var,1,clr1)

    # plotting the embedded cells
    for cell in embedded_cell_pos
        cell = equidistant_normal_points_matrix(cell, 1.0)
        CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=4, color=:black)
        CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=3, color=:red)
    end

    # plotting the final interface
    plotInterface_no_cb!(ax, u,var,length(u),clr2)

end

function plotGrowthSimulation_no_cb_emb_cell_orientation!(ax, u, var, embedded_cell_pos, clr1, clr2)

    #println(size(embedded_cell_pos))
    orientations = KubaPhD.calc_cell_orientation(embedded_cell_pos)[3]
    #println(length(orientations))
    #println(length(embedded_cell_pos))

    # plotting the initial interface
    plotInterface_no_cb!(ax,u,var,1, clr1)

    # plotting the embedded cells
    for ii in eachindex(embedded_cell_pos)
        cell = equidistant_normal_points_matrix(embedded_cell_pos[ii], 1.0)
        println(orientations[ii])
        CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=3, color=:black)
        CairoMakie.lines!(ax, cell[1,:], cell[2,:], linewidth=2, color=ones(length(cell[2,:]))*orientations[ii], colormap=:Spectral, colorrange=(0,2π))
    end

    # plotting the final interface
    plotInterface_no_cb!(ax, u,var,length(u),clr2)

end


function plotCellCount!(ax_cell, ax_embedded, t, cell_count, embed_cell_count)
    CairoMakie.lines!(ax_cell, t, cell_count, color=:blue, linewidth=2)
    CairoMakie.lines!(ax_embedded, t, embed_cell_count, color=:red, linewidth=2)
end

function plotCellCount!(ax_cell, ax_embedded, t, cell_count, embed_cell_count, clr1, clr2)
    CairoMakie.lines!(ax_cell, t, cell_count, color=clr1, linewidth=2)
    CairoMakie.lines!(ax_embedded, t, embed_cell_count, color=clr2, linewidth=2)
end

function plotInterfaceDensity!(ax, t, den)
    mean_interface_density = zeros(length(t))
    for i in 1:length(t)
        mean_interface_density[i] = sum(den[i]) / length(den[i])
    end
    CairoMakie.lines!(ax, t, den[1][1].*ones(size(t)), color=:black, linestyle=:dash, linewidth=2)
    CairoMakie.scatter!(ax, t, mean_interface_density, color=:blue, markersize=5)
    CairoMakie.limits!(ax, 0, t[end]+0.1, 0.0, ceil(maximum(mean_interface_density),digits=1))
end

function plotInterfaceDensity!(ax, t, den, clr)
    mean_interface_density = zeros(length(t))
    for i in 1:length(t)
        mean_interface_density[i] = sum(den[i]) / length(den[i])
    end
    CairoMakie.lines!(ax, t, den[1][1].*ones(size(t)), color=:black, linestyle=:dash, linewidth=2)
    CairoMakie.scatter!(ax, t, mean_interface_density, color=clr, markersize=5)
    CairoMakie.limits!(ax, 0, t[end]+0.1, 0.0, ceil(maximum(mean_interface_density),digits=1))
end

function plotTissueCoverage!(ax, t, area)
    CairoMakie.lines!(ax, t, area ./ area[1], color=:black, linewidth=2)
end

function plotTissueCoverage!(ax, t, area, clr)
    CairoMakie.lines!(ax, t, area ./ area[1], color=clr, linewidth=3)
end

function plotTissueArea!(ax, t, area)
    CairoMakie.lines!(ax, t, area[1] .- area, color=:black, linewidth=2)
end

function plotTissueArea!(ax, t, area, clr)
    CairoMakie.lines!(ax, t, area[1] .- area, color=clr, linewidth=2)
end

function plotEmbedmentRate!(ax, t, Ot, embed_cell_rates, kf)
    CairoMakie.lines!(ax, t, ones(length(t)).*Ot, color=:black, linestyle=:dash, linewidth=2)
    CairoMakie.lines!(ax, t, embed_cell_rates, color=:red, linewidth=2)
end

function plotEmbedDensity!(ax, t, embed_count, Ω, CellMech, Embed)
    # calculate embed density = #/Ω
    embed_density = embed_count ./ (Ω[1] .- Ω)
    CairoMakie.lines!(ax, t, (ones(size(t)).*Embed.rate) ./ CellMech.kf, linewidth=5, linestyle=:dash, color=:black)
    CairoMakie.lines!(ax, t, embed_density, linewidth = 5, color=:red)
end

function plotCellOrientation!(ax, θ, R, clr)
    CairoMakie.stairs!(ax, [θ; θ[1]], [R; R[1]],linewidth=3,color=clr)
end

function plotCellOrientation_at_t!(ax, all_embedded_cells, embedded_cell_count, idx, clr)
    θ,R = calc_cell_orientation(all_embedded_cells[end - Int64((embedded_cell_count[1][idx])[1]) + 1:end])
    plotCellOrientation!(ax, θ, R, clr)
end

function plotEmbedCellDensityMap!(ax, sol)

    init_interface = sol.u[1]
    final_interface = sol.u[end]

    init_interface2 = zeros(size(init_interface,2)+1, size(init_interface,1))
    init_interface2[1:end-1,1] =  init_interface[1,:]
    init_interface2[1:end-1,2] =  init_interface[2,:]
    init_interface2[end,1] = init_interface[1,1]; init_interface2[end,2] = init_interface[2,1]; 

    final_interface2 = zeros(size(final_interface,2)+1, size(final_interface,1))
    final_interface2[1:end-1,1] =  final_interface[1,:]
    final_interface2[1:end-1,2] =  final_interface[2,:]
    final_interface2[end,1] = final_interface[1,1]; final_interface2[end,2] = final_interface[2,1]; 
    reverse!(final_interface2, dims=1)

    boundary_points = [
        [[(p[1], p[2]) for p in eachrow(init_interface2)]],
        [[(p[1], p[2]) for p in eachrow(final_interface2)]],
    ]

    boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
    tri = triangulate(points; boundary_nodes)
    refine!(tri; max_area=1e-3get_area(tri))

    x_coords = Float64[]
    y_coords = Float64[]
    for point in tri.points
        push!(x_coords, point[1])
        push!(y_coords, point[2])
    end

    z = zeros(size(x_coords,1))

    cells = KubaPhD.MCP(embedded_cells);

    for cell in cells
        try
            z[DelaunayTriangulation.get_nearest_neighbour(tri,cell)] = z[DelaunayTriangulation.get_nearest_neighbour(tri,cell)] + 1.0
        catch err
            continue
        end
    end

    tr = tricontourf!(ax, tri, z, colormap = Reverse(:dense), levels=30)
    return tr
end

function plotContinuumSimulation!(ax, θ, R, ρ, cmap, cbar_limits)
    for i in 1:5000:size(R,1)
        CairoMakie.lines!(ax,[R[i,:]; R[i,1]].*cos.([θ;θ[1]]), [R[i,:]; R[i,1]].*sin.([θ;θ[1]]) ,color=[ρ[i,:];ρ[i,1]],linewidth=4, colorrange=cbar_limits, colormap=cmap)
    end
end

function plotMeanInterfaceDensity_DiscVsCont!(ax, all_disc_solutions, cont_sol)
    # discrete simulations 
    all_mean_interface_density = Vector{Float64}[];
    for sim in all_disc_solutions
        mean_interface_density = zeros(size(sim.t,1))
        for time in eachindex(sim.t)
            mean_interface_density[time] = sum(sim.Density[time])/size(sim.Density[time],1)
        end
        push!(all_mean_interface_density, mean_interface_density)
    end
    # columns are different simulations, rows are time points
    all_mean_interface_density_matrix = hcat(all_mean_interface_density...);

    total_mean = sum(all_mean_interface_density_matrix,dims=2)./size(all_mean_interface_density_matrix,2)

    # continuum simulations
    mean_int_density_cont = zeros(size(cont_sol.ρ,1))
    for time in eachindex(mean_int_density_cont)
        mean_int_density_cont[time] = sum(cont_sol.ρ[time])/size(cont_sol.ρ[time],1)
    end

    CairoMakie.lines!(ax,cont_sol.t, mean_int_density_cont ,color=:midnightblue, linestyle=:solid,linewidth=4, label="Continuum")
    CairoMakie.lines!(ax,all_disc_solutions[1].t, [total_mean...] ,color=:springgreen3, linestyle=:dash,linewidth=4, label="100 realisation average")
    labels = ["Realisation 1", "Realisation 2", "Realisation 3", "Realisation 4", "Realisation 5"]
    indecies = [1,26,51,76,100]
    ii = 1
    for col in eachcol(all_mean_interface_density_matrix[:,1:5])
       CairoMakie.scatter!(ax, all_disc_solutions[1].t[indecies], col[indecies], markersize=10, strokewidth = 2, strokecolor = :black, label = labels[ii])
       ii += 1
    end

end

function plotMeanInterfaceDensity_DiscVsCont!(ax1, ax2, all_disc_solutions_1, all_disc_solutions_2, all_disc_solutions_3, cont_sol)
    # discrete simulations 
    all_mean_interface_density_list = []
    total_mean_list = []
    for all_disc_solutions in [all_disc_solutions_1, all_disc_solutions_2, all_disc_solutions_3]
        all_mean_interface_density = Vector{Float64}[];
        for sim in all_disc_solutions
            mean_interface_density = zeros(size(sim.t,1))
            for time in eachindex(sim.t)
                mean_interface_density[time] = sum(sim.Density[time])/size(sim.Density[time],1)
            end
            push!(all_mean_interface_density, mean_interface_density)
        end
        # columns are different simulations, rows are time points
        all_mean_interface_density_matrix = hcat(all_mean_interface_density...);
        total_mean = sum(all_mean_interface_density_matrix,dims=2)./size(all_mean_interface_density_matrix,2)
        
        push!(all_mean_interface_density_list, all_mean_interface_density_matrix)
        push!(total_mean_list, total_mean)
    end

    all_mean_interface_density_matrix_1 = all_mean_interface_density_list[1]
    total_mean_1 = total_mean_list[1]
    total_mean_2 = total_mean_list[2]
    total_mean_3 = total_mean_list[3]

    # continuum simulations
    mean_int_density_cont = zeros(size(cont_sol.ρ,1))
    for time in eachindex(mean_int_density_cont)
        mean_int_density_cont[time] = sum(cont_sol.ρ[time])/size(cont_sol.ρ[time],1)
    end

    CairoMakie.lines!(ax1,cont_sol.t, mean_int_density_cont ,color=:midnightblue, linestyle=:solid,linewidth=4, label="Continuum")
    CairoMakie.lines!(ax1,all_disc_solutions_1[1].t, [total_mean_1...] ,color=:black, linestyle=:dash,linewidth=4, label="m = 2 average")
    CairoMakie.lines!(ax1,all_disc_solutions_2[1].t, [total_mean_2...] ,color=:darkorange, linestyle=:dash,linewidth=4, label="m = 4 average")
    CairoMakie.lines!(ax1,all_disc_solutions_2[1].t, [total_mean_3...] ,color=:springgreen3, linestyle=:dash,linewidth=4, label="m = 10 average")
    labels = ["Realisation 1", "Realisation 2", "Realisation 3", "Realisation 4", "Realisation 5"]
    indecies = [1,26,51,76,100]
    ii = 1
    for col in eachcol(all_mean_interface_density_matrix_1[:,1:5])
       CairoMakie.scatter!(ax1, all_disc_solutions_1[1].t[indecies], col[indecies], markersize=10, strokewidth = 2, strokecolor = :black) #, label = labels[ii])
       ii += 1
    end


end

function plotTissueArea_DiscVsCont!(ax, all_disc_solutions, cont_sol)
    # discrete simulations 
    all_tissue_area = Vector{Float64}[];
    for sim in all_disc_solutions
        push!(all_tissue_area, sim.Ω[1] .- sim.Ω)
    end
    # columns are different simulations, rows are time points
    all_tissue_area_matrix = hcat(all_tissue_area...);

    total_area = sum(all_tissue_area_matrix,dims=2)./size(all_tissue_area_matrix,2)

    # continuum simulations
    total_area_cont = cont_sol.R[:,1] .^2 .* π

    CairoMakie.lines!(ax,cont_sol.t, total_area_cont[1] .- total_area_cont ,color=:midnightblue, linestyle=:solid,linewidth=4, label="Continuum")
    CairoMakie.lines!(ax,all_disc_solutions[1].t, [total_area...] ,color=:springgreen3, linestyle=:dash,linewidth=4, label="Discrete average")
    labels = ["Realisation 1", "Realisation 2", "Realisation 3", "Realisation 4", "Realisation 5"]
    indecies = [1,26,51,76,100]
    ii = 1
    for col in eachcol(all_tissue_area_matrix[:,1:5])
       CairoMakie.scatter!(ax, all_disc_solutions[1].t[indecies], col[indecies], markersize=10, strokewidth = 2, strokecolor = :black, label = labels[ii])
       ii += 1
    end
end

function plotTissueArea_DiscVsCont!(ax, all_disc_solutions_1, all_disc_solutions_2, all_disc_solutions_3, cont_sol)
    all_tissue_area_list = []
    total_area_list = []
    for all_disc_solutions in [all_disc_solutions_1, all_disc_solutions_2, all_disc_solutions_3]
        all_tissue_area = Vector{Float64}[];
        for sim in all_disc_solutions
            push!(all_tissue_area, sim.Ω[1] .- sim.Ω)
        end
        # columns are different simulations, rows are time points
        all_tissue_area_matrix = hcat(all_tissue_area...);
        total_area = sum(all_tissue_area_matrix,dims=2)./size(all_tissue_area_matrix,2)
        
        push!(all_tissue_area_list, all_tissue_area_matrix)
        push!(total_area_list, total_area)
    end

    all_tissue_area_matrix_1 = all_tissue_area_list[1]
    total_area_1 = total_area_list[1]
    total_area_2 = total_area_list[2]
    total_area_3 = total_area_list[3]

    # continuum simulations
    total_area_cont = cont_sol.R[:,1] .^2 .* π

    CairoMakie.lines!(ax,cont_sol.t, total_area_cont[1] .- total_area_cont ,color=:midnightblue, linestyle=:solid,linewidth=4, label="Continuum")
    CairoMakie.lines!(ax,all_disc_solutions_1[1].t, [total_area_1...] ,color=:black, linestyle=:dash,linewidth=4, label="m = 2 average")
    CairoMakie.lines!(ax,all_disc_solutions_2[1].t, [total_area_2...] ,color=:darkorange, linestyle=:dash,linewidth=4, label="m = 4 average")
    CairoMakie.lines!(ax,all_disc_solutions_3[1].t, [total_area_3...] ,color=:springgreen3, linestyle=:dash,linewidth=4, label="m = 10 average")
    labels = ["Realisation 1", "Realisation 2", "Realisation 3", "Realisation 4", "Realisation 5"]
    indecies = [1,26,51,76,100]
    ii = 1
    for col in eachcol(all_tissue_area_matrix_1[:,1:5])
       CairoMakie.scatter!(ax, all_disc_solutions_1[1].t[indecies], col[indecies], markersize=10, strokewidth = 2, strokecolor = :black, label = labels[ii])
       ii += 1
    end
end