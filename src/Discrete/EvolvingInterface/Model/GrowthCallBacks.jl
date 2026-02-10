# Cell behaviour Callback

"""
    event_affect!(integrator)

Update the state of `integrator` based on probabilistic cellular events.

This function modifies `integrator` in place. It uses the parameters and state from `integrator` to compute probabilities for different cellular events: proliferation (`prolif`), death (`death`), and embedding (`embed`). Based on these probabilities and random draws, it updates the state of `integrator`.

# Arguments
- `integrator`: The integrator object containing the current state and parameters. The parameters are expected to be a tuple containing:

# Returns
`nothing`. The function modifies `integrator` in place.
"""
function event_affect!(integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = integrator.p
    u = integrator.u
    curr_t = integrator.t
    (p, a, e, pe) = cell_probs(u, curr_t, Domain.m, SimTime, Prolif, Death, Embed, ProlifEmbed)
    (r1, r2, r3) = rand(3)
    total_prob = sum(p) + sum(a) + sum(e) + sum(pe)
    
    if r1 < total_prob # check if event has occurred
        event_type = r2 * total_prob
        event_occuring = ""
        if event_type < sum(p) # prolif occurs
            # finding proliferating cell index
            idx = find_cell_index(p, r3 * sum(p))
            left_spring_index = idx * Domain.m - (Domain.m - 1)
            right_spring_index = idx * Domain.m + 1
            right_spring_index = right_spring_index > size(u, 2) ? 1 : right_spring_index
            # interpolating & inserting new daughter cell spring boundaries
            centre_spring_pos = Domain.m % 2 == 0 ? u[:, Int64(left_spring_index + Domain.m / 2)] : (u[:, Int64(left_spring_index + (Domain.m - 1) / 2)] + u[:, Int64(left_spring_index + (Domain.m - 1) / 2 + 1)]) / 2
            interp_points_left = LinearInterp(u[:, left_spring_index], centre_spring_pos, Domain.m)
            interp_points_right = LinearInterp(centre_spring_pos, u[:, right_spring_index], Domain.m)
            for i = 1:Domain.m - 1
                deleteat!(u, left_spring_index + 1)
            end
            new_prolif_spring_pos = vcat(interp_points_left, [centre_spring_pos], interp_points_right)
            for point in reverse(new_prolif_spring_pos)
                insert!(u, left_spring_index + 1, point)
            end
            # inserting new daughter cell spring mechanical properties
            for i in 1:Domain.m
                insert!(CellMech.kₛ, left_spring_index + 1, CellMech.kₛ[idx])
                insert!(CellMech.kf, left_spring_index + 1, CellMech.kf[idx])
                insert!(CellMech.a, left_spring_index + 1, CellMech.a[idx])
                insert!(CellMech.growth_dir, left_spring_index + 1, CellMech.growth_dir[idx])
            end
            event_occuring = "prolif"
        elseif event_type < sum(p) + sum(a) # death occurs
            idx = find_cell_index(a, r3 * sum(a))
            spring_index = idx * Domain.m - (Domain.m - 1)
            boundaries_midpoint = Domain.m % 2 == 0 ? u[:, Int64(spring_index + Domain.m / 2)] : (u[:, Int64(spring_index + (Domain.m - 1) / 2)] + u[:, Int64(spring_index + (Domain.m - 1) / 2 + 1)]) / 2
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

        elseif event_type < sum(p) + sum(a) + sum(e) # embed only occurs
            #println("embed occured @ t = ", integrator.t)
            idx = find_cell_index(e, r3 * sum(e))
            spring_index = idx * Domain.m - (Domain.m - 1)
            store_embedded_cell(u, spring_index, Domain.m, integrator.t)
            boundaries_midpoint = Domain.m % 2 == 0 ? u[:, Int64(spring_index + Domain.m / 2)] : (u[:, Int64(spring_index + (Domain.m - 1) / 2)] + u[:, Int64(spring_index + (Domain.m - 1) / 2 + 1)]) / 2
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
            #resize!(integrator, (2, size(integrator.u, 2)))
        else # prolif and embed simultaneously occur 
            idx = find_cell_index(pe, r3 * sum(pe))
            spring_index = idx * Domain.m - (Domain.m - 1)
            store_embedded_cell(u, spring_index, Domain.m, integrator.t)
        end
        resize!(integrator, (2, size(integrator.u, 2)))

        if size(integrator.u, 2) != size(CellMech.kₛ, 1) 
            error("sizing error with: ", event_occuring)
        end
    end
    nothing
end

function event_affect_old!(integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = integrator.p
    u = integrator.u
    curr_t = integrator.t
    (p, a, e, pe) = cell_probs(u, curr_t, Domain.m, SimTime, Prolif, Death, Embed, ProlifEmbed)
    (r1, r2, r3) = rand(3)
    total_prob = sum(p) + sum(a) + sum(e) + sum(pe)
    
    if r1 < total_prob # check if event has occurred
        event_type = r2 * total_prob
        event_occuring = ""
        if event_type < sum(p) # prolif occurs
            # finding proliferating cell index
            idx = find_cell_index(p, r3 * sum(p))
            left_spring_index = idx * Domain.m - (Domain.m - 1)
            right_spring_index = idx * Domain.m + 1
            right_spring_index = right_spring_index > size(u, 2) ? 1 : right_spring_index
            # interpolating & inserting new daughter cell spring boundaries
            centre_spring_pos = Domain.m % 2 == 0 ? u[:, Int64(left_spring_index + Domain.m / 2)] : (u[:, Int64(left_spring_index + (Domain.m - 1) / 2)] + u[:, Int64(left_spring_index + (Domain.m - 1) / 2 + 1)]) / 2
            interp_points_left = LinearInterp(u[:, left_spring_index], centre_spring_pos, Domain.m)
            interp_points_right = LinearInterp(centre_spring_pos, u[:, right_spring_index], Domain.m)
            for i = 1:Domain.m - 1
                deleteat!(u, left_spring_index + 1)
            end
            new_prolif_spring_pos = vcat(interp_points_left, [centre_spring_pos], interp_points_right)
            for point in reverse(new_prolif_spring_pos)
                insert!(u, left_spring_index + 1, point)
            end
            # inserting new daughter cell spring mechanical properties
            for i in 1:Domain.m
                insert!(CellMech.kₛ, left_spring_index + 1, CellMech.kₛ[idx])
                insert!(CellMech.kf, left_spring_index + 1, CellMech.kf[idx])
                insert!(CellMech.a, left_spring_index + 1, CellMech.a[idx])
                insert!(CellMech.growth_dir, left_spring_index + 1, CellMech.growth_dir[idx])
            end
            event_occuring = "prolif"
        elseif event_type < sum(p) + sum(a) # death occurs
            idx = find_cell_index(a, r3 * sum(a))
            spring_index = idx * Domain.m - (Domain.m - 1)
            #boundaries_midpoint = Domain.m % 2 == 0 ? u[:, Int64(spring_index + Domain.m / 2)] : (u[:, Int64(spring_index + (Domain.m - 1) / 2)] + u[:, Int64(spring_index + (Domain.m - 1) / 2 + 1)]) / 2
            for i = 1:Domain.m
                deleteat!(u, spring_index)
                deleteat!(CellMech.kₛ, spring_index)
                deleteat!(CellMech.kf, spring_index)
                deleteat!(CellMech.a, spring_index)
                deleteat!(CellMech.growth_dir, spring_index)
            end
            #if spring_index > size(u, 2)
            #    u[:, 1] .= boundaries_midpoint
            #else
            #    u[:, spring_index] .= boundaries_midpoint
            #end
            event_occuring = "death"

        elseif event_type < sum(p) + sum(a) + sum(e) # embed only occurs
            #println("embed occured @ t = ", integrator.t)
            idx = find_cell_index(e, r3 * sum(e))
            spring_index = idx * Domain.m - (Domain.m - 1)
            store_embedded_cell(u, spring_index, Domain.m, integrator.t)
            #boundaries_midpoint = Domain.m % 2 == 0 ? u[:, Int64(spring_index + Domain.m / 2)] : (u[:, Int64(spring_index + (Domain.m - 1) / 2)] + u[:, Int64(spring_index + (Domain.m - 1) / 2 + 1)]) / 2
            for i = 1:Domain.m
                deleteat!(u, spring_index)
                deleteat!(CellMech.kₛ, spring_index)
                deleteat!(CellMech.kf, spring_index)
                deleteat!(CellMech.a, spring_index)
                deleteat!(CellMech.growth_dir, spring_index)
            end
            #if spring_index > size(u, 2)
            #    u[:, 1] .= boundaries_midpoint
            #else
            #    u[:, spring_index] .= boundaries_midpoint
            #end
            event_occuring = "embed"
            #resize!(integrator, (2, size(integrator.u, 2)))
        else # prolif and embed simultaneously occur 
            idx = find_cell_index(pe, r3 * sum(pe))
            spring_index = idx * Domain.m - (Domain.m - 1)
            store_embedded_cell(u, spring_index, Domain.m, integrator.t)
        end
        resize!(integrator, (2, size(integrator.u, 2)))

        if size(integrator.u, 2) != size(CellMech.kₛ, 1) 
            error("sizing error with: ", event_occuring)
        end
    end
    nothing
end

function LinearInterp(a::ElasticVector, b::ElasticVector, m::Int64)
    dt = 1 / m
    T = dt:dt:1-dt
    Interp_P = [a + t * (b - a) for t in T]
    return Interp_P
end



terminate_affect!(integrator) = terminate!(integrator)

function density_lim_condition(u,t,integrator)
    (m,kₛ,η,kf,l₀,δt,growth_dir,domain_type,btype,prolif,death,embed,α,β,γ,q_lim) = integrator.p
    q = calc_cell_densities(u,m)
    flag = sum(q .> q_lim)
    return flag > 0.0 ? true : false
end

# Callback for stopping at a certain area coverage
function area_lim_condition(u,t,integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = integrator.p
    curr_tissue_area_coverage = (Domain.Ω_0 - Ω(u)) / Domain.Ω_0;
    if curr_tissue_area_coverage > 0.90
        println("Simulation terminated at t = ", integrator.t)
    end
    return curr_tissue_area_coverage > 0.90
end

# Callback for stopping if interface overlaps
function interface_overlap_condition(u,t,integrator)
    if interface_overlaps(integrator.u)
        println("Interface overlaps. Simulation terminated at t = ", integrator.t)
    end
    return interface_overlaps(integrator.u)
end

# supporting functions
function segments_intersect(a1, a2, b1, b2)
    cross(v1, v2) = v1[1]*v2[2] - v1[2]*v2[1]
    d1 = a2 - a1
    d2 = b2 - b1
    delta = b1 - a1

    denom = cross(d1, d2)
    if denom == 0
        return false
    end

    t = cross(delta, d2) / denom
    u = cross(delta, d1) / denom

    return 0 < t < 1 && 0 < u < 1
end

function interface_overlaps(points)
    N = size(points, 2)
    for i in 1:N
        a1 = points[:, i]
        a2 = points[:, mod1(i + 1, N)]  # mod1 handles wrapping (1-based indexing)

        for j in i+1:N
            # Skip adjacent or identical (wraparound) edges
            if abs(i - j) <= 1 || (i == 1 && j == N)
                continue
            end
            b1 = points[:, j]
            b2 = points[:, mod1(j + 1, N)]
            if segments_intersect(a1, a2, b1, b2)
                return true
            end
        end
    end
    return false
end