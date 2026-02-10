"""
    δ(rᵢ₊₁, rᵢ)

Calculate the Euclidean distance between two points `rᵢ₊₁` and `rᵢ`.

# Arguments
- `rᵢ₊₁`: The first point in space.
- `rᵢ`: The second point in space.

# Returns
The Euclidean distance between the two points.
"""
#δ(rᵢ₊₁, rᵢ) = .√(sum((rᵢ₊₁ - rᵢ).^2,dims=size(rᵢ)))

function δ(rᵢ₊₁, rᵢ)
    val = zeros(1,size(rᵢ,2))
    val .= .√(sum((rᵢ₊₁- rᵢ).^2,dims=1))
    return val
end

"""
    ρ(rᵢ₊₁, rᵢ)

Calculate the reciprocal of the distance (interpreted as density) between two points `rᵢ₊₁` and `rᵢ`.

# Arguments
- `rᵢ₊₁`: The first point in space.
- `rᵢ`: The second point in space.

# Returns
The reciprocal of the distance between the two points.
"""
ρ(rᵢ₊₁, rᵢ) = 1 ./ δ(rᵢ₊₁, rᵢ);

"""
    τ(rᵢ₊₁, rᵢ₋₁)

Calculate the unit tangent vector between two neighboring points `rᵢ₊₁` and `rᵢ₋₁`.

# Arguments
- `rᵢ₊₁`: The point after the central point in space.
- `rᵢ₋₁`: The point before the central point in space.

# Returns
The unit tangent vector between the two points.
"""
τ(rᵢ₊₁, rᵢ₋₁) = (rᵢ₊₁ - rᵢ₋₁) ./ δ(rᵢ₊₁, rᵢ₋₁)

"""
    n(rᵢ₊₁, rᵢ₋₁, type)

Calculate the unit normal vector at a point `rᵢ` between two neighboring points `rᵢ₊₁` and `rᵢ₋₁`. The orientation of the normal vector depends on the specified `type`.

# Arguments
- `rᵢ₊₁`: The point after the central point in space.
- `rᵢ₋₁`: The point before the central point in space.
- `type`: A string specifying the orientation of the normal vector, either "inward" or any other value for outward orientation.

# Returns
The unit normal vector at the point.
"""
function n(rᵢ₊₁, rᵢ₋₁,type) 
    if type == "inward"
        -oftype(τ(rᵢ₊₁, rᵢ₋₁),vcat(transpose.([-τ(rᵢ₊₁, rᵢ₋₁)[:,2], τ(rᵢ₊₁, rᵢ₋₁)[:,1]])...)')
    else
        oftype(τ(rᵢ₊₁, rᵢ₋₁),vcat(transpose.([-τ(rᵢ₊₁, rᵢ₋₁)[:,2], τ(rᵢ₊₁, rᵢ₋₁)[:,1]])...)')
    end
end


function calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir)
    l = .√(sum((rᵢ₊₁ - rᵢ).^2, dims=1)) # calculating length
    ρ = 1.0 ./ l # calculating density
    dv = (rᵢ₊₁ - rᵢ) ./ l # directional vector
    τ = (rᵢ₊₁ - rᵢ₋₁) ./ .√(sum((rᵢ₊₁ - rᵢ₋₁).^2, dims=1)) # tangent vector

    if growth_dir == "inward" # normal vector
        n = [-1.0, 1.0] .* circshift(dv, 1)
    else
        n = [1.0, -1.0] .* circshift(dv, 1)
    end

    return l, ρ, dv, τ, n # returning all calculated values
end

function calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir::ElasticVector{String})
    l = .√(sum((rᵢ₊₁ - rᵢ).^2, dims=1)) # calculating length
    ρ = 1.0 ./ l # calculating density
    dv = (rᵢ₊₁ - rᵢ) ./ l # directional vector
    τ = (rᵢ₊₁ - rᵢ₋₁) ./ .√(sum((rᵢ₊₁ - rᵢ₋₁).^2, dims=1)) # tangent vector

    n = zeros(size(dv))
    for i in eachindex(growth_dir)
        shift_dv = circshift(dv[:, i], 1)
        if growth_dir[i] == "inward"
            n[:, i] = [-1.0, 1.0] .* shift_dv
        else
            n[:, i] = [1.0, -1.0] .* shift_dv
        end
    end

    return l, ρ, dv, τ, n # returning all calculated values
end

"""
    lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ)

Implementation of geometric solution for moving a cell in the normal direction. Calculate the intersection points of two lines defined by points `rₘ₁`, `rₗ` and `rₘ₂`, `rᵣ`.

This function computes the intersection point of each pair of lines. If the lines are parallel (the determinant is zero), it returns the midpoint of `rₘ₁` and `rₘ₂`. Otherwise, it calculates the intersection point using the parameters `u` and `t`. If the intersection lies within the segments, the intersection point is returned; otherwise, the midpoint is returned.

# Arguments
- `rₘ₁`: Starting point of the first line segment.
- `rₗ`: Ending point of the first line segment.
- `rₘ₂`: Starting point of the second line segment.
- `rᵣ`: Ending point of the second line segment.

# Returns
A vector of intersection points for each pair of line segments.
"""
function lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ)
    intersect = zeros(size(rₘ₁))
    r = rₗ .- rₘ₁
    s = rᵣ .- rₘ₂
    d = r[1, :] .* s[2, :] .- r[2, :] .* s[1, :]

    for i in eachindex(d)
        if d[i] == 0
            @views intersect[:, i] .= (rₘ₁[:, i] .+ rₘ₂[:, i]) ./ 2
        else
            u = ((rₘ₂[1, i] - rₘ₁[1, i]) * r[2, i] - (rₘ₂[2, i] - rₘ₁[2, i]) * r[1, i]) / d[i]
            t = ((rₘ₂[1, i] - rₘ₁[1, i]) * s[2, i] - (rₘ₂[2, i] - rₘ₁[2, i]) * s[1, i]) / d[i]

            if 0 ≤ u ≤ 1 && 0 ≤ t ≤ 1
                @views intersect[:, i] .= rₘ₁[:, i] .+ t .* r[:, i]
            else
                @views intersect[:, i] .= (rₘ₁[:, i] .+ rₘ₂[:, i]) ./ 2
            end
        end
    end
    return intersect
end

"""
    Ω(p)

Calculate the area of a polygon defined by points in `p`.

This function computes the area using the shoelace formula. The polygon is defined by a set of points `p`, and the function iterates through these points to calculate the area.

# Arguments
- `p`: A matrix where each row represents a point of the polygon in 2D space.

# Returns
The absolute area of the polygon.
"""
function Ω(p)
    A = 0
    for ii in axes(p,2)
        if ii == size(p,2)
            A += (p[1,ii]*p[2,1] -  p[2,ii]*p[1,1])
        else
            A += (p[1,ii]*p[2,ii+1] -  p[2,ii]*p[1,ii+1])
        end
    end
    return abs(A)/2;
end

"""
    ωκ(rᵢ₋₁, rᵢ, rᵢ₊₁)

Calculate the areas of triangles formed by consecutive triplets of points.

This function takes three points `rᵢ₋₁`, `rᵢ`, and `rᵢ₊₁` and calculates the area of the triangle formed by these points. It is typically used in a loop over a sequence of points to calculate the area of each triangle formed by consecutive triplets.

# Arguments
- `rᵢ₋₁`: The first point of the triangle.
- `rᵢ`: The second point of the triangle.
- `rᵢ₊₁`: The third point of the triangle.

# Returns
A vector of areas for each triangle formed by the input points.
"""
function ωκ(rᵢ₋₁, rᵢ, rᵢ₊₁)
    triVector = [rᵢ₋₁' rᵢ' rᵢ₊₁']
    A = zeros(size(triVector,1))
    for ii in axes(triVector,1)
        A[ii] = Ω(reshape(triVector[ii,:],(2,3))')
    end
    return A
end


"""
    κ(rᵢ₋₁, rᵢ, rᵢ₊₁)

Approximate the curvature of a shape using the Menger method. This method is based on the areas of triangles formed by consecutive triplets of points.

# Arguments
- `rᵢ₋₁`: The point before the current point in space.
- `rᵢ`: The current point in space.
- `rᵢ₊₁`: The point after the current point in space.

# Returns
The approximated curvature at the point `rᵢ`.

# Reference
Anoshkina, Elena V., Alexander G. Belyaev, and Hans-Peter Seidel. "Asymtotic Analysis of Three-Point Approximations of Vertex Normals and Curvatures." VMV. 2002.
"""
function κ(rᵢ₋₁, rᵢ, rᵢ₊₁)

    A = ωκ(rᵢ₋₁,rᵢ,rᵢ₊₁)
    l1 = .√(sum((rᵢ- rᵢ₋₁).^2,dims=1))
    l2 = .√(sum((rᵢ₊₁- rᵢ).^2,dims=1))
    l3 = .√(sum((rᵢ₊₁- rᵢ₋₁).^2,dims=1))

    return (4*A)./(l1.*l2.*l3)'
end

function approx_cell_κ(cell)
    mean_κ = 0
    sub_κ = size(cell,2)-2

    if sub_κ == 0
        return Inf
    else
        for ii = 1:sub_κ
            mean_κ = mean_κ .+ κ(cell[:,ii]', cell[:,ii+1]', cell[:,ii+2]')
        end
        mean_κ = mean_κ/sub_κ
    end
    return mean_κ
end

function calc_shape_centroid(data, idx)
    # Calculate the centroid of the shape
    x = data.u[idx][1, :]
    y = data.u[idx][2, :]
    A = data.Ω[idx]

    x_centroid = 1 / (6 * A) * sum((x[1:end-1] + x[2:end]) .* (x[1:end-1] .* y[2:end] - x[2:end] .* y[1:end-1]))
    y_centroid = 1 / (6 * A) * sum((y[1:end-1] + y[2:end]) .* (x[1:end-1] .* y[2:end] - x[2:end] .* y[1:end-1]))

    return [x_centroid, y_centroid]
end

