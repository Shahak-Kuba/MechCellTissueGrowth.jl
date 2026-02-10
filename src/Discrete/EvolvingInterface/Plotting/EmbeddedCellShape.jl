

function equidistant_normal_points_matrix(segments::Matrix{Float64}, distance::Float64)
    # find the total length of the cell
    
    points_above = []
    points_below = []

    # Loop through each segment (each column in the matrix)
    for i in 1:size(segments, 2)-1
        p1 = (segments[1, i], segments[2, i])         # Starting point
        p2 = (segments[1, i+1], segments[2, i+1])     # Ending point

        # Check if the segment is valid
        if i == size(segments, 2)
            break  # Skip the last point as it has no next point
        end

        # Calculate the midpoint
        midpoint = ((p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2)

        # Calculate the direction vector of the line segment
        direction = (p2[1] - p1[1], p2[2] - p1[2])

        # Calculate the length of the direction vector
        length = sqrt(direction[1]^2 + direction[2]^2)

        # Normalize the direction vector
        if length == 0
            throw(ArgumentError("The points p1 and p2 must not be the same"))
        end
        normalized_direction = (direction[1] / length, direction[2] / length)

        # Calculate the normal vector (perpendicular to the direction)
        normal = (-normalized_direction[2], normalized_direction[1])

        # Calculate the two points above and below the midpoint
        # density dependent size
        #point_above = [midpoint[1] + normal[1] * (1/length) * distance, midpoint[2] + normal[2] * (1/length) * distance]
        #point_below = [midpoint[1] - normal[1] * (1/length) * distance, midpoint[2] - normal[2] * (1/length) * distance]
        # constant size
        point_above = [midpoint[1] + normal[1] * distance, midpoint[2] + normal[2] * distance]
        point_below = [midpoint[1] - normal[1] * distance, midpoint[2] - normal[2] * distance]

        # Add points to the array
        push!(points_above, point_above)
        push!(points_below, point_below)
    end

    ordered_points = [hcat(segments[:,1], hcat(points_above...), segments[:,end], hcat(reverse(points_below)...), segments[:,1])][1]

    return ordered_points
end

"""
# Example usage
segments = embedded_cells[1]

distance = 10.0

p = Plots.plot()
for cell in embedded_cells
    points_array = equidistant_normal_points_matrix(cell, distance)
    p = Plots.plot!(points_array[1,:],points_array[2,:],linewidth=7, legend=false, color=:black)
    p = Plots.plot!(points_array[1,:],points_array[2,:],linewidth=2, legend=false, color=:red)
end

display(p)

points_array = equidistant_normal_points_matrix(segments, distance)

Plots.plot(segments[1,:],segments[2,:],linewidth=10)
Plots.scatter!(segments[1,:], segments[2,:])
Plots.scatter!(points_array[1,:], points_array[2,:])
Plots.plot!(points_array[1,:],points_array[2,:],linewidth=10)
"""


function equidistant_normal_points_matrix(segments::Matrix{Float64}, distance::Float64, cell_length_thresh::Float64)
    # find the total length of the cell
    segment_lengths = Float64[]
    for i in 1:size(segments, 2)-1
        p1 = (segments[1, i], segments[2, i])         #
        p2 = (segments[1, i+1], segments[2, i+1])     # Ending point
        segment_length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
        push!(segment_lengths, segment_length)
    end
    cs = cumsum(segment_lengths)
    cell_body_segments = count(<(cell_length_thresh), cs)
    #println("cell length: $cs")
    
    if cell_body_segments >= size(segments, 2)-1 # if the cell is too short, we make it longer for visualisation
        scaled_pts, s, Lold = scale_curve_to_length(segments, cs[end], cell_length_thresh; anchor=:centroid)
        segments = scaled_pts
        println("scaled cell from length $Lold to $cell_length_thresh with scale factor $s")
    end

    points_above = []
    points_below = []

    shape = :ellipse
    println("cell body shape: $shape")

    if shape == :ellipse
        KubaPhD.build_outline_ellipse!(points_above, points_below, segments; max_radius=distance)
    else
        # Loop through each segment (each column in the matrix)
        for i in 1:min(cell_body_segments, size(segments, 2)-1)
            p1 = (segments[1, i], segments[2, i])         # Starting point
            p2 = (segments[1, i+1], segments[2, i+1])     # Ending point

            # Check if the segment is valid
            if i == size(segments, 2)
                break  # Skip the last point as it has no next point
            end

            # Calculate the midpoint
            midpoint = ((p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2)

            # Calculate the direction vector of the line segment
            direction = (p2[1] - p1[1], p2[2] - p1[2])

            # Calculate the length of the direction vector
            length = sqrt(direction[1]^2 + direction[2]^2)

            # Normalize the direction vector
            if length == 0
                throw(ArgumentError("The points p1 and p2 must not be the same"))
            end
            normalized_direction = (direction[1] / length, direction[2] / length)

            # Calculate the normal vector (perpendicular to the direction)
            normal = (-normalized_direction[2], normalized_direction[1])

            # Calculate the two points above and below the midpoint
            # density dependent size
            #point_above = [midpoint[1] + normal[1] * (1/length) * distance, midpoint[2] + normal[2] * (1/length) * distance]
            #point_below = [midpoint[1] - normal[1] * (1/length) * distance, midpoint[2] - normal[2] * (1/length) * distance]
            # constant size
            point_above = [midpoint[1] + normal[1] * distance, midpoint[2] + normal[2] * distance]
            point_below = [midpoint[1] - normal[1] * distance, midpoint[2] - normal[2] * distance]

            # Add points to the array
            push!(points_above, point_above)
            push!(points_below, point_below)
        end
    end

    ordered_points = [hcat(segments[:,1], hcat(points_above...), segments[:,cell_body_segments+1], hcat(reverse(points_below)...), segments[:,1])][1]

    return ordered_points
end



"""
Compute polyline length: sum of Euclidean distances between consecutive points.
pts is an N×2 (or N×d) matrix, each row a point.

function polyline_length(pts::AbstractMatrix)
    N = size(pts, 1)
    @assert N ≥ 2
    L = 0.0
    for i in 1:(N-1)
        L += norm(view(pts, i+1, :) .- view(pts, i, :))
    end
    return L
end

"""

"""
Scale a curve to have target length Lnew, preserving shape.
anchor = :centroid, :first, :last, or an explicit vector.
Returns (scaled_pts, scale_factor, old_length).
"""
function scale_curve_to_length(pts::AbstractMatrix, Lold, Lnew::Real; anchor=:centroid)
    @assert Lold > 0
    s = Lnew / Lold

    a = if anchor === :centroid
        vec(mean(pts, dims=2))
    elseif anchor === :first
        vec(pts[:, 1])
    elseif anchor === :last
        vec(pts[end, :])
    elseif anchor isa AbstractVector
        anchor
    else
        error("anchor must be :centroid, :first, :last, or a vector")
    end

    scaled = similar(pts, float(eltype(pts)))
    for i in 1:size(pts,2)
        scaled[:, i] = a .+ s .* (vec(pts[:, i]) .- a)
    end
    return scaled, s, Lold
end


# segments: 2×N matrix, columns are points
# max_radius: half-width at the fattest point
function build_outline_ellipse!(points_above, points_below, segments; max_radius)
    n = size(segments, 2)
    n < 2 && return

    # Arc-length parameterization so taper is smooth even with uneven spacing
    s = zeros(Float64, n)
    for i in 2:n
        dx = segments[1,i] - segments[1,i-1]
        dy = segments[2,i] - segments[2,i-1]
        s[i] = s[i-1] + hypot(dx, dy)
    end
    total = s[end]
    total == 0 && throw(ArgumentError("All segment points are identical"))

    for i in 1:(n-1)
        p1x, p1y = segments[1,i],   segments[2,i]
        p2x, p2y = segments[1,i+1], segments[2,i+1]

        dx, dy = p2x - p1x, p2y - p1y
        len = hypot(dx, dy)
        len == 0 && continue

        tx, ty = dx/len, dy/len                  # tangent (unit)
        nx, ny = -ty, tx                         # normal (unit)
        mx, my = (p1x+p2x)/2, (p1y+p2y)/2         # midpoint

        # normalized position along the body (0 at start, 1 at end), using midpoint's s
        u = ((s[i] + s[i+1]) / 2) / total

        # Ellipse cross-section profile: r(u) = b * sqrt(1 - (2u-1)^2)
        # peaks at u=0.5, goes to 0 at ends
        x = 2u - 1
        r = max_radius * sqrt(max(0.0, 1 - x^2))

        push!(points_above, [mx + nx*r, my + ny*r])
        push!(points_below, [mx - nx*r, my - ny*r])
    end
end