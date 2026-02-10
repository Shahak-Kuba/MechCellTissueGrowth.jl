

function interpolate_polygon_matrix(V::Matrix{Float64}, M::Int)
    # Ensure it's 2 x N
    @assert size(V, 1) == 2 "Input matrix must be 2 x N"

    N = size(V, 2)
    closed = all(V[:, 1] .== V[:, end])
    if !closed
        V = hcat(V, V[:, 1])  # close the loop
        N += 1
    end

    # Compute edge lengths
    segment_lengths = [norm(V[:, i+1] - V[:, i]) for i in 1:N-1]
    total_length = sum(segment_lengths)
    spacing = total_length / M

    points = Matrix{Float64}(undef, 2, M)
    current_pos = 0.0
    segment_index = 1
    seg_len = segment_lengths[1]

    for i in 1:M
        target_dist = (i - 1) * spacing

        # Advance to the correct segment
        while current_pos + seg_len < target_dist && segment_index < N - 1
            current_pos += seg_len
            segment_index += 1
            seg_len = segment_lengths[segment_index]
        end

        t = (target_dist - current_pos) / seg_len
        p1 = V[:, segment_index]
        p2 = V[:, segment_index + 1]
        points[:, i] = (1 - t) * p1 + t * p2
    end

    return points  # 2 x M matrix of equidistant points
end

function find_closest_point_in_centroid_direction(sol_0::Matrix{Float64}, sol_end::Matrix{Float64}, centroid::Vector{Float64})
    closest_pairs = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, size(sol_0, 2))
    for i in 1:size(sol_0, 2)
        p1 = sol_0[:, i]
        line_vec = centroid - p1
        distances = Vector{Tuple{Float64, Vector{Float64}}}()
        
        for j in 1:size(sol_end, 2)
            p2 = sol_end[:, j]
            # Project p2 onto the line formed by p1 and centroid
            proj = p1 + (dot(p2 - p1, line_vec) / dot(line_vec, line_vec)) * line_vec
            # Compute the perpendicular distance from p2 to the line
            dist_to_line = norm(p2 - proj)
            push!(distances, (dist_to_line, p2))
        end
        
        # Sort distances and take the closest 10 points
        sorted_distances = sort(distances, by=x -> x[1])
        closest_10 = sorted_distances[1:min(10, length(sorted_distances))]
        
        # Find the point with the shortest Euclidean distance to p1
        min_dist = Inf
        closest_point = nothing
        for (_, p2) in closest_10
            dist = norm(p2 - p1)
            if dist < min_dist
                min_dist = dist
                closest_point = p2
            end
        end
        
        if closest_point === nothing
            closest_pairs[i] = (p1, [NaN, NaN])  # Assign a default value if no closest point is found
        else
            closest_pairs[i] = (p1, closest_point)
        end
    end
    return closest_pairs
end

function atan2pi(y, x)
    angle = atan(y, x)
    return angle >= 0 ? angle : angle + 2π
end

function closest_points_to_line(interface_0, interface_end, p1, p2)
    sol_u = interface_0
    output = interface_end
    # Ensure the input is valid
    @assert size(sol_u, 1) == 2 "sol_u must be a 2 x N matrix"
    @assert size(output, 1) == 2 "output must be a 2 x M matrix"

    # Line direction vector
    line_dir = normalize(p2 - p1)

    # Compute the angle of the line with the x-axis
    theta = atan2pi(line_dir[2], line_dir[1])

    # Function to compute the perpendicular distance from a point to the line
    function point_to_line_distance(point)
        projection = p1 + dot(point - p1, line_dir) * line_dir
        return norm(point - projection)
    end

    # Function to check if a point makes the desired angle with the x-axis
    function is_in_dir_of_line(point1, point2, theta)
        line_dir = normalize(point1 - point2)
        theta1 = atan2pi(line_dir[2], line_dir[1])
        return abs(theta1 - theta) < 0.1  # Allow a small tolerance
    end

    # Find the 10 closest points in sol_u to the line that satisfy the angle condition
    distances_sol_u = [(point_to_line_distance(sol_u[:, i]), i) for i in 1:size(sol_u, 2)]
    closest_2_sol_u = sort(distances_sol_u, by=x -> x[1])[1:min(20, length(distances_sol_u))]

    # Find the 10 closest points in output to the line that satisfy the angle condition
    distances_output = [(point_to_line_distance(output[:, i]), i) for i in 1:size(output, 2)]
    closest_2_output = sort(distances_output, by=x -> x[1])[1:min(20, length(distances_output))]

    # Find the pair of points (one from each interface) with the shortest distance between them
    min_dist = Inf
    closest_point_sol_u = nothing
    closest_point_output = nothing

    for (_, idx_u) in closest_2_sol_u
        for (_, idx_out) in closest_2_output
            dist = norm(sol_u[:, idx_u] - output[:, idx_out])
            if dist < min_dist && is_in_dir_of_line(sol_u[:, idx_u], output[:, idx_out], theta)
                min_dist = dist
                closest_point_sol_u = sol_u[:, idx_u]
                closest_point_output = output[:, idx_out]
                #println("Closest point in sol_u: ", idx_u)
                #println("Closest point in output: ", idx_out)
            end
        end
    end

    return min_dist, closest_point_sol_u, closest_point_output
end

function find_largest_WTh_asymmetry_ratio(interface_0, interface_end, centroid)
    points_in_der_centr= find_closest_point_in_centroid_direction(interface_0::Matrix{Float64}, interface_end::Matrix{Float64}, centroid::Vector{Float64})
    WTh_asymmetry_ratio = 1.0
    p1_o = nothing
    p2_o = nothing
    p3_o = nothing
    p4_o = nothing
    for points in points_in_der_centr
        p1 = points[1]
        p2 = points[2]
        # Calculate the distance between p1 and p2
        min_dist = norm(p1 - p2)
        # Calculate the angle between the line segment and the x-axis
        max_dist, p3, p4 = closest_points_to_line(interface_0, interface_end, p1, p2)
        # Calculate the WTh asymmetry ratio
        if min_dist / max_dist < WTh_asymmetry_ratio
            p1_o = p1
            p2_o = p2
            p3_o = p3
            p4_o = p4
            WTh_asymmetry_ratio = min_dist / max_dist
        end
    end

    return WTh_asymmetry_ratio, p1_o, p2_o, p3_o, p4_o
end

function Generate_rectangles_in_growth_region(p1, p2, p3, p4, width)

    function generate_rectangle(p1, p2, width)
        midpoint = (p1 + p2) / 2
        direction = normalize(p2 - p1)
        perpendicular = [-direction[2], direction[1]]  # Rotate 90 degrees
        half_width = width / 2

        # Compute rectangle corners
        corner1 = midpoint - half_width * perpendicular + (p1 - midpoint)
        corner2 = midpoint + half_width * perpendicular + (p1 - midpoint)
        corner3 = midpoint + half_width * perpendicular + (p2 - midpoint)
        corner4 = midpoint - half_width * perpendicular + (p2 - midpoint)

        # Draw the rectangle
        #CairoMakie.poly!(ax, [corner1[1], corner2[1], corner3[1], corner4[1], corner1[1]],
                            #[corner1[2], corner2[2], corner3[2], corner4[2], corner1[2]],
                            #color=color, strokewidth=0)

        # Calculate the area of the rectangle
        length = norm(p2 - p1)
        area = length * width

        return (corner1, corner2, corner3, corner4), area
    end

    corners_small, area_small = generate_rectangle(p1, p2, width)
    corners_large, area_large = generate_rectangle(p3, p4, width)

    return (corners_small, area_small), (corners_large, area_large)

end

function count_embedded_cells_in_rectangles(embedded_cells, small_rect, large_rect)
    function is_point_in_rectangle(point, corners)
        # Use a point-in-polygon check for rectangles
        x, y = point
        x_coords = [corner[1] for corner in corners]
        y_coords = [corner[2] for corner in corners]
        n = length(corners)
        inside = false
        j = n
        for i in 1:n
            xi, yi = x_coords[i], y_coords[i]
            xj, yj = x_coords[j], y_coords[j]
            if (yi > y) != (yj > y) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
                inside = !inside
            end
            j = i
        end
        return inside
    end

    small_corners = small_rect[1]
    large_corners = large_rect[1]

    count_small = 0
    count_large = 0

    for cell in embedded_cells
        # Check if any point in the cell is within the small rectangle
        if any(point -> is_point_in_rectangle(point, small_corners), eachcol(cell))
            count_small += 1
        end
        # Check if any point in the cell is within the large rectangle
        if any(point -> is_point_in_rectangle(point, large_corners), eachcol(cell))
            count_large += 1
        end
    end

    return count_small, count_large
end

function analyse_wallThickness_density(all_solutions, all_embedded_cell_pos)
    results = DataFrame(
        idx = Int[],
        small_area = Float64[],
        large_area = Float64[],
        count_small = Int[],
        count_large = Int[],
        ξ_small = Float64[],
        ξ_large = Float64[],
        WThRatio = Float64[]
    )

    for idx in 1:length(all_solutions)
        sol = all_solutions[idx]
        embedded_cells = all_embedded_cell_pos[idx]
        sol_0 = interpolate_polygon_matrix(sol.u[1], 2000)
        sol_end = interpolate_polygon_matrix(sol.u[end], 2000)

        centroid = KubaPhD.calc_shape_centroid(sol, size(sol.u, 1))
        WTh_asymmetry_ratio, p1, p2, p3, p4 = find_largest_WTh_asymmetry_ratio(sol_0, sol_end, centroid)


        small_rect, large_rect = Generate_rectangles_in_growth_region(p1, p2, p3, p4, 40)
        count_small, count_large = count_embedded_cells_in_rectangles(embedded_cells, small_rect, large_rect)

        ξ_small = count_small / small_rect[2]
        ξ_large = count_large / large_rect[2]

        push!(results, (
            idx = idx,
            small_area = small_rect[2],
            large_area = large_rect[2],
            count_small = count_small,
            count_large = count_large,
            ξ_small = ξ_small,
            ξ_large = ξ_large,
            WThRatio = WTh_asymmetry_ratio
        ))
    end

    return results
end