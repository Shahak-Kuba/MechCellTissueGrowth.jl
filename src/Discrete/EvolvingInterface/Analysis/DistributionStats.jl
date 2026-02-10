"""
Mean Cell Position
"""
function MCP(points)
    return sum.(points,dims=2)./size.(points,2)
end

"""
Point Correlation

for reference see https://physics.emory.edu/faculty/weeks/idl/gofr2.html
"""
function PairCorrelation(points,r,dr,Λ)
    g = Float64[]
    # Get the mean positions of all embedded cells
    spatial_points = MCP(points);
    # looping over all radius values
    for rad in r
        N = 0;
        # looping over all points
        for point in spatial_points
            p_rep = repeat(point,1,size(spatial_points,1))
            # measure the distance between all other cells and the specific cell
            all_distances = KubaPhD.δ(p_rep', hcat(spatial_points...)');
            N = N + sum((rad .< all_distances) .&& (all_distances .< (rad+dr)))
        end
        push!(g, (N/size(spatial_points,1))/(2π*rad*dr))
    end
    return g
end

"""
Ripley's K function

"""

function RKfunc(points, r)
    k = Float64[]
    # Get the mean positions of all embedded cells
    spatial_points = MCP(points);
    # looping over all radius values
    for rad in r
        N = 0;
        # looping over all points
        for point in spatial_points
            p_rep = repeat(point,1,size(spatial_points,1))
            # measure the distance between all other cells and the specific cell
            all_distances = KubaPhD.δ(p_rep', hcat(spatial_points...)');
            N = N + sum((rad .>= all_distances))
        end
        push!(k, (2π*rad^2))/(size(spatial_points,1)*(size(spatial_points,1)-1) * N)
    end
    return k
end

"""
Mean K-Nearest Neighbour

"""

function meanKNN(points)
    mKNN = 0.0;
    kNearNeighbour = Float64[];
    spatial_points = MCP(points);
    for point in spatial_points
        p_rep = repeat(point,1,size(spatial_points,1))
        # measure the distance between all other cells and the specific cell
        all_distances = KubaPhD.δ(p_rep', hcat(spatial_points...)');
        push!(kNearNeighbour, minimum(all_distances[all_distances .> 0.0]))
    end
    mKNN = sum(kNearNeighbour)/size(spatial_points,1)
    return mKNN
end

"""
Radial Cell Distribution

"""

function RadDensity(points, r, dr)
    k = Float64[]
    # Get the mean positions of all embedded cells
    spatial_points = MCP(points)
    # looping over all radius values
    for rad in r
        N = 0
        # measure the distance between the origin and all cells
        all_distances = KubaPhD.δ([0.0, 0.0], hcat(spatial_points...)')
        N = sum((rad .< all_distances) .& (all_distances .< (rad + dr)))
        push!(k, N / (2π * rad * dr))
    end
    return k
end

"""
    Wall thickness density
    
"""


function WallDensity(embedded_count, Ω)
    Ω_upp = Ω
    Ω_low = circshift(Ω,-1)
    Ω_diff = (Ω_low - Ω_upp)[1:end-1]

    N_upp = embedded_count
    N_low = circshift(embedded_count,-1)
    N_diff = (N_low - N_upp)[1:end-1]

    cellDensity = N_diff./Ω_diff
    #while size(cellDensity,1) < max_size
    #    push!(cellDensity, -1000)
    #    push!(Ω_diff, -1000)
    #    push!(N_diff, -1000)
    #end
    return cellDensity, Ω_diff, N_diff
end





