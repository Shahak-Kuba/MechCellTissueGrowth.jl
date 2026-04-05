# Cell behaviour Callbacks

"""
    event_affect!(integrator)

Update the state of `integrator` based on probabilistic cellular events.
Uses `SimulationState` from `p` instead of global variables.
Resizes `ODECache` after cell count changes.
"""
function event_affect!(integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, cache, state = integrator.p
    u = integrator.u
    curr_t = integrator.t
    (p_prob, a_prob, e_prob, pe_prob) = cell_probs(u, curr_t, Domain.m, SimTime, Prolif, Death, Embed, ProlifEmbed)
    (r1, r2, r3) = rand(3)
    total_prob = sum(p_prob) + sum(a_prob) + sum(e_prob) + sum(pe_prob)

    if r1 < total_prob
        event_type = r2 * total_prob
        event_occuring = ""

        if event_type < sum(p_prob)
            # ── Proliferation ──
            idx = find_cell_index(p_prob, r3 * sum(p_prob))
            left_spring_index = idx * Domain.m - (Domain.m - 1)
            right_spring_index = idx * Domain.m + 1
            right_spring_index = right_spring_index > size(u, 2) ? 1 : right_spring_index

            centre_spring_pos = Domain.m % 2 == 0 ?
                u[:, Int64(left_spring_index + Domain.m / 2)] :
                (u[:, Int64(left_spring_index + (Domain.m - 1) / 2)] +
                 u[:, Int64(left_spring_index + (Domain.m - 1) / 2 + 1)]) / 2

            interp_points_left  = LinearInterp(u[:, left_spring_index], centre_spring_pos, Domain.m)
            interp_points_right = LinearInterp(centre_spring_pos, u[:, right_spring_index], Domain.m)

            for i = 1:Domain.m - 1
                deleteat!(u, left_spring_index + 1)
            end
            new_prolif_spring_pos = vcat(interp_points_left, [centre_spring_pos], interp_points_right)
            for point in reverse(new_prolif_spring_pos)
                insert!(u, left_spring_index + 1, point)
            end

            # Insert new daughter cell spring mechanical properties
            for i in 1:Domain.m
                insert!(CellMech.kₛ, left_spring_index + 1, CellMech.kₛ[idx])
                insert!(CellMech.kf, left_spring_index + 1, CellMech.kf[idx])
                insert!(CellMech.a, left_spring_index + 1, CellMech.a[idx])
                insert!(CellMech.growth_dir, left_spring_index + 1, CellMech.growth_dir[idx])
            end
            event_occuring = "prolif"

        elseif event_type < sum(p_prob) + sum(a_prob)
            # ── Death ──
            idx = find_cell_index(a_prob, r3 * sum(a_prob))
            spring_index = idx * Domain.m - (Domain.m - 1)
            boundaries_midpoint = Domain.m % 2 == 0 ?
                u[:, Int64(spring_index + Domain.m / 2)] :
                (u[:, Int64(spring_index + (Domain.m - 1) / 2)] +
                 u[:, Int64(spring_index + (Domain.m - 1) / 2 + 1)]) / 2

            for i = 1:Domain.m
                deleteat!(u, spring_index)
                deleteat!(CellMech.kₛ, spring_index)
                deleteat!(CellMech.kf, spring_index)
                deleteat!(CellMech.a, spring_index)
                deleteat!(CellMech.growth_dir, spring_index)
            end
            if spring_index > size(u, 2)
                u[:, 1] .= boundaries_midpoint
            else
                u[:, spring_index] .= boundaries_midpoint
            end
            event_occuring = "death"

        elseif event_type < sum(p_prob) + sum(a_prob) + sum(e_prob)
            # ── Embedding ──
            idx = find_cell_index(e_prob, r3 * sum(e_prob))
            spring_index = idx * Domain.m - (Domain.m - 1)
            store_embedded_cell(state, u, spring_index, Domain.m, integrator.t)
            boundaries_midpoint = Domain.m % 2 == 0 ?
                u[:, Int64(spring_index + Domain.m / 2)] :
                (u[:, Int64(spring_index + (Domain.m - 1) / 2)] +
                 u[:, Int64(spring_index + (Domain.m - 1) / 2 + 1)]) / 2

            for i = 1:Domain.m
                deleteat!(u, spring_index)
                deleteat!(CellMech.kₛ, spring_index)
                deleteat!(CellMech.kf, spring_index)
                deleteat!(CellMech.a, spring_index)
                deleteat!(CellMech.growth_dir, spring_index)
            end
            if spring_index > size(u, 2)
                u[:, 1] .= boundaries_midpoint
            else
                u[:, spring_index] .= boundaries_midpoint
            end
            event_occuring = "embed"

        else
            # ── Proliferation + Embedding simultaneously ──
            idx = find_cell_index(pe_prob, r3 * sum(pe_prob))
            spring_index = idx * Domain.m - (Domain.m - 1)
            store_embedded_cell(state, u, spring_index, Domain.m, integrator.t)
        end

        resize!(integrator, (2, size(integrator.u, 2)))
        # Resize the ODE cache to match the new number of nodes
        resize_cache!(cache, size(integrator.u, 2))

        if size(integrator.u, 2) != size(CellMech.kₛ, 1)
            error("sizing error with: ", event_occuring)
        end
    end
    nothing
end


function LinearInterp(a::ElasticVector, b::ElasticVector, m::Int64)
    dt = 1.0 / m
    T  = dt:dt:1-dt
    return [a + t * (b - a) for t in T]
end


terminate_affect!(integrator) = terminate!(integrator)

function density_lim_condition(u, t, integrator)
    Domain = integrator.p[1]
    CellMech = integrator.p[2]
    q = calc_cell_densities(u, Domain.m)
    q_lim = 100.0  # default limit; adjust as needed
    return sum(q .> q_lim) > 0
end

# Callback for stopping at a certain area coverage
function area_lim_condition(u, t, integrator)
    Domain = integrator.p[1]
    curr_tissue_area_coverage = (Domain.Ω_0 - Ω(u)) / Domain.Ω_0
    if curr_tissue_area_coverage > 0.90
        println("Simulation terminated at t = ", integrator.t)
    end
    return curr_tissue_area_coverage > 0.90
end

# Callback for stopping if interface overlaps
function interface_overlap_condition(u, t, integrator)
    if interface_overlaps(integrator.u)
        println("Interface overlaps. Simulation terminated at t = ", integrator.t)
    end
    return interface_overlaps(integrator.u)
end

# Supporting functions
function segments_intersect(a1, a2, b1, b2)
    cross(v1, v2) = v1[1] * v2[2] - v1[2] * v2[1]
    d1 = a2 - a1
    d2 = b2 - b1
    delta = b1 - a1

    denom = cross(d1, d2)
    denom == 0 && return false

    t = cross(delta, d2) / denom
    u = cross(delta, d1) / denom

    return 0 < t < 1 && 0 < u < 1
end

function interface_overlaps(points)
    N = size(points, 2)
    for i in 1:N
        a1 = points[:, i]
        a2 = points[:, mod1(i + 1, N)]
        for j in i+1:N
            (abs(i - j) <= 1 || (i == 1 && j == N)) && continue
            b1 = points[:, j]
            b2 = points[:, mod1(j + 1, N)]
            segments_intersect(a1, a2, b1, b2) && return true
        end
    end
    return false
end
