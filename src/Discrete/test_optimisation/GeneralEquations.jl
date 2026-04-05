
# =====================================================================
#  Scalar / inline helpers for the hot path (zero-allocation)
# =====================================================================

"""Scalar force computation — used inside the ODE loop."""
@inline function force_scalar(l::Float64, kₛ::Float64, a::Float64, ft::ForceType)
    if ft == Hookean
        return kₛ * (l - a)
    elseif ft == Nonlinear
        return kₛ * (1.0 / a - 1.0 / l)
    else  # Nonlinear2
        return kₛ * a * a * (1.0 / a - 1.0 / l)
    end
end

"""Scalar line-intersection for two 2D line segments — used inside the ODE loop."""
@inline function _line_intersect(
    rm1x::Float64, rm1y::Float64, rlx::Float64, rly::Float64,
    rm2x::Float64, rm2y::Float64, rrx::Float64, rry::Float64
)
    rx = rlx - rm1x
    ry = rly - rm1y
    sx = rrx - rm2x
    sy = rry - rm2y
    d  = rx * sy - ry * sx

    if d == 0.0
        return (rm1x + rm2x) * 0.5, (rm1y + rm2y) * 0.5
    end

    dmx = rm2x - rm1x
    dmy = rm2y - rm1y
    u_p = (dmx * ry - dmy * rx) / d
    t_p = (dmx * sy - dmy * sx) / d

    if 0.0 ≤ u_p ≤ 1.0 && 0.0 ≤ t_p ≤ 1.0
        return rm1x + t_p * rx, rm1y + t_p * ry
    else
        return (rm1x + rm2x) * 0.5, (rm1y + rm2y) * 0.5
    end
end

"""
    colwise_dot!(result, A, B)

Compute `result[j] = dot(A[:, j], B[:, j])` for each column.
Replaces the costly `diag(A' * B)` pattern.
"""
function colwise_dot!(result, A, B)
    @inbounds for j in axes(A, 2)
        s = 0.0
        for i in axes(A, 1)
            s += A[i, j] * B[i, j]
        end
        result[j] = s
    end
    nothing
end

"""Zero-allocation cumulative-sum search (replaces `find_cell_index`)."""
function find_cell_index(arr::AbstractVector{Float64}, threshold::Float64)
    s = 0.0
    @inbounds for i in eachindex(arr)
        s += arr[i]
        s >= threshold && return i
    end
    return length(arr) - 1
end

# =====================================================================
#  Vectorised helpers (used by post-processing & analysis code)
# =====================================================================

"""Euclidean distance between columns of two 2×N matrices, returned as 1×N."""
function δ(rᵢ₊₁, rᵢ)
    return .√(sum((rᵢ₊₁ .- rᵢ).^2, dims=1))
end

"""Density (1/distance) between columns of two 2×N matrices."""
ρ(rᵢ₊₁, rᵢ) = 1.0 ./ δ(rᵢ₊₁, rᵢ)

"""Unit tangent vector between two neighbouring points."""
τ(rᵢ₊₁, rᵢ₋₁) = (rᵢ₊₁ .- rᵢ₋₁) ./ δ(rᵢ₊₁, rᵢ₋₁)

"""Unit normal vector (2D rotation of tangent)."""
function n(rᵢ₊₁, rᵢ₋₁, type)
    t = τ(rᵢ₊₁, rᵢ₋₁)
    if type == "inward" || type == GD_INWARD
        return vcat(-t[2:2, :], t[1:1, :])
    else
        return vcat(t[2:2, :], -t[1:1, :])
    end
end

"""
    calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir)

Compute spring lengths, densities, directional vectors, tangent vectors, and normals.
Dispatches on scalar `GrowthDir` or per-spring `AbstractVector{Int8}`.
"""
function calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir::GrowthDir)
    l  = .√(sum((rᵢ₊₁ .- rᵢ).^2, dims=1))
    ρ  = 1.0 ./ l
    dv = (rᵢ₊₁ .- rᵢ) ./ l
    τ  = (rᵢ₊₁ .- rᵢ₋₁) ./ .√(sum((rᵢ₊₁ .- rᵢ₋₁).^2, dims=1))
    gs = Int8(growth_dir)
    n  = similar(dv)
    @inbounds for j in axes(dv, 2)
        n[1, j] =  gs * dv[2, j]
        n[2, j] = -gs * dv[1, j]
    end
    return l, ρ, dv, τ, n
end

function calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir::AbstractVector{Int8})
    l  = .√(sum((rᵢ₊₁ .- rᵢ).^2, dims=1))
    ρ  = 1.0 ./ l
    dv = (rᵢ₊₁ .- rᵢ) ./ l
    τ  = (rᵢ₊₁ .- rᵢ₋₁) ./ .√(sum((rᵢ₊₁ .- rᵢ₋₁).^2, dims=1))
    n  = similar(dv)
    @inbounds for j in eachindex(growth_dir)
        gs = growth_dir[j]
        n[1, j] =  gs * dv[2, j]
        n[2, j] = -gs * dv[1, j]
    end
    return l, ρ, dv, τ, n
end

# Backward-compat: accept a plain string for FB code
function calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir::String)
    gd = growth_dir == "inward" ? GD_INWARD : GD_OUTWARD
    return calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, gd)
end

"""
    lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ)

Vectorised line-intersection for 2×N point sets.
Used by post-processing; the ODE hot path uses `_line_intersect` instead.
"""
function lineIntersection(rₘ₁, rₗ, rₘ₂, rᵣ)
    intersect = similar(rₘ₁)
    r = rₗ .- rₘ₁
    s = rᵣ .- rₘ₂
    d = r[1, :] .* s[2, :] .- r[2, :] .* s[1, :]

    @inbounds for i in eachindex(d)
        if d[i] == 0
            intersect[1, i] = (rₘ₁[1, i] + rₘ₂[1, i]) / 2
            intersect[2, i] = (rₘ₁[2, i] + rₘ₂[2, i]) / 2
        else
            u = ((rₘ₂[1, i] - rₘ₁[1, i]) * r[2, i] - (rₘ₂[2, i] - rₘ₁[2, i]) * r[1, i]) / d[i]
            t = ((rₘ₂[1, i] - rₘ₁[1, i]) * s[2, i] - (rₘ₂[2, i] - rₘ₁[2, i]) * s[1, i]) / d[i]
            if 0 ≤ u ≤ 1 && 0 ≤ t ≤ 1
                intersect[1, i] = rₘ₁[1, i] + t * r[1, i]
                intersect[2, i] = rₘ₁[2, i] + t * r[2, i]
            else
                intersect[1, i] = (rₘ₁[1, i] + rₘ₂[1, i]) / 2
                intersect[2, i] = (rₘ₁[2, i] + rₘ₂[2, i]) / 2
            end
        end
    end
    return intersect
end

"""Shoelace-formula area of polygon with columns as points."""
function Ω(p)
    A = 0.0
    N = size(p, 2)
    @inbounds for ii in 1:N
        jj = ii == N ? 1 : ii + 1
        A += p[1, ii] * p[2, jj] - p[2, ii] * p[1, jj]
    end
    return abs(A) / 2
end

"""Triangle areas from consecutive triplets of points."""
function ωκ(rᵢ₋₁, rᵢ, rᵢ₊₁)
    triVector = [rᵢ₋₁' rᵢ' rᵢ₊₁']
    A = zeros(size(triVector, 1))
    for ii in axes(triVector, 1)
        A[ii] = Ω(reshape(triVector[ii, :], (2, 3))')
    end
    return A
end

"""Menger curvature approximation."""
function κ(rᵢ₋₁, rᵢ, rᵢ₊₁)
    A  = ωκ(rᵢ₋₁, rᵢ, rᵢ₊₁)
    l1 = .√(sum((rᵢ .- rᵢ₋₁).^2, dims=1))
    l2 = .√(sum((rᵢ₊₁ .- rᵢ).^2, dims=1))
    l3 = .√(sum((rᵢ₊₁ .- rᵢ₋₁).^2, dims=1))
    return (4 .* A) ./ (l1 .* l2 .* l3)'
end

function approx_cell_κ(cell)
    sub_κ = size(cell, 2) - 2
    sub_κ == 0 && return Inf
    mean_κ = 0.0
    for ii in 1:sub_κ
        mean_κ += κ(cell[:, ii]', cell[:, ii+1]', cell[:, ii+2]')[1]
    end
    return mean_κ / sub_κ
end

function calc_shape_centroid(data, idx)
    x = data.u[idx][1, :]
    y = data.u[idx][2, :]
    A = data.Ω[idx]
    cross_terms = x[1:end-1] .* y[2:end] .- x[2:end] .* y[1:end-1]
    x_centroid = sum((x[1:end-1] .+ x[2:end]) .* cross_terms) / (6A)
    y_centroid = sum((y[1:end-1] .+ y[2:end]) .* cross_terms) / (6A)
    return [x_centroid, y_centroid]
end
