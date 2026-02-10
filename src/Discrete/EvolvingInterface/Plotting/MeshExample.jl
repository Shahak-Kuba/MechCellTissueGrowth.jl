# Run simulation first

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

fig = Figure(fontsize = 42, size = (800, 800))
wh = (width = 600, height = 600)
ax = Axis(fig[1, 1]; wh...)
triplot!(ax, tri)
#hidedecorations!(ax)
#hidespines!(ax)
fig

save("artwork2.pdf",fig)

x_coords = Float64[]
y_coords = Float64[]
for point in tri.points
    push!(x_coords, point[1])
    push!(y_coords, point[2])
end

z = zeros(size(x_coords,1))

cells = KubaPhD.MCP(embed_pos);


for cell in cells
    try
        z[DelaunayTriangulation.get_nearest_neighbour(tri,cell)] = z[DelaunayTriangulation.get_nearest_neighbour(tri,cell)] + 1.0
    catch err
        continue
    end
end


fig = Figure(fontsize = 42, size = (800, 800))
wh = (width = 600, height = 600)
ax = Axis(fig[1, 1]; wh...)

tr = tricontourf!(ax, tri, z, colormap = Reverse(:dense), levels=30;)

fig, ax, tr = tricontourf(tri, z, colormap = Reverse(:dense), levels=30; 
        axis = (; aspect = 1, title = "Embedded cell heatmap", titlefont = "Times New Roman", titlesize=20, xlabelfont="Times New Roman", xlabelsize=20, ylabelfont="Times New Roman", ylabelsize=20, xlabel = "x [μm]", ylabel = "y [μm]"))
Colorbar(fig[1, 2], tr; width = 20, label = "Number of embedded cells", labelsize=20, labelfont = "Times New Roman", labelpadding = 10.0)
#hidedecorations!(ax)
#hidespines!(ax)
fig

