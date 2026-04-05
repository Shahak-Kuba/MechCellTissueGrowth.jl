"""
    PostCalcs1D(u, p, CellMech_at_t)

Post-calculation for 1D simulation data.
"""
function PostCalcs1D(u, p, CellMech_at_t)
    Domain  = p[1]
    CellMech = p[2]
    SimTime = p[3]

    kₘ         = CellMech_at_t[1]
    aₘ         = CellMech_at_t[2]
    kfₘ        = CellMech_at_t[3]
    growth_dir = CellMech_at_t[4]
    η          = CellMech.η

    ∑F      = zeros(size(u, 1))
    density = zeros(size(u, 1))
    vₙ      = zeros(size(u, 1))
    ψ       = zeros(size(u, 1))
    Κ       = zeros(size(u, 1))

    density .= 1.0 ./ (Domain.m .* diff(u[1, :]))

    return ∑F, vₙ, density, ψ, Κ
end


"""
    PostCalcs2D(u, p, CellMech_at_t)

Post-calculation for 2D simulation data.
Uses `colwise_dot!` instead of `diag(A' * B)`.
Computes Vₙ only once (fixes double-computation bug).
"""
function PostCalcs2D(u, p, CellMech_at_t)
    Domain  = p[1]
    CellMech = p[2]
    SimTime = p[3]

    kₘ         = CellMech_at_t[1]
    aₘ         = CellMech_at_t[2]
    kfₘ        = CellMech_at_t[3]
    growth_dir = CellMech_at_t[4]
    η          = CellMech.η

    N = size(u, 2)
    ∑F      = zeros(N)
    density = zeros(N)
    vₙ      = zeros(N)
    ψ       = zeros(N)
    Κ       = zeros(N)

    uᵢ₊₁ = circshift(u, (0, -1))
    uᵢ₋₁ = circshift(u, (0, 1))

    # Convert growth_dir to Int8 vector if needed for dispatch
    gd = growth_dir isa AbstractVector{Int8} ? growth_dir : Int8.(growth_dir)
    l, ρ, dv, τ, n_vec = calc_l_ρ_dv_τ_n(uᵢ₊₁, u, uᵢ₋₁, gd)

    # Force computation: use colwise_dot! instead of diag(A' * B)
    force_field = Fₛ⁺(l, kₘ, aₘ, CellMech.restoring_force) .* dv .+
                  Fₛ⁻(circshift(l, (0, 1)), kₘ, aₘ, CellMech.restoring_force) .* circshift(dv, (0, 1))
    colwise_dot!(∑F, force_field, τ)

    density .= ((circshift(ρ, (0, -1)) .+ ρ) ./ (2 * Domain.m))'
    ψ .= ∑F ./ (kₘ .* aₘ)
    Κ .= κ(uᵢ₋₁, u, uᵢ₊₁)

    # Compute Vₙ once and extract x/y components
    Vn_full = Vₙ(uᵢ₋₁, u, uᵢ₊₁, ρ, n_vec, kfₘ, SimTime.δt)
    vₙx = @view Vn_full[1, :]
    vₙy = @view Vn_full[2, :]
    vₙ .= .√(vₙx.^2 .+ vₙy.^2)

    return ∑F, vₙ, density, ψ, Κ
end


"""
    postSimulation(sol, p, AllCellMech)

Process simulation solution into `SimResults_t`.
"""
function postSimulation(sol, p, AllCellMech)
    Domain  = p[1]
    SimTime = p[3]

    c = size(sol.t, 1)
    s = min(c, size(AllCellMech, 1))
    Area       = Vector{Float64}(undef, s)
    Cell_Count = Vector{Float64}(undef, s)
    ∑F      = Vector{Vector{Float64}}(undef, 0)
    ψ       = Vector{Vector{Float64}}(undef, 0)
    DENSITY = Vector{Vector{Float64}}(undef, 0)
    vₙ_out  = Vector{Vector{Float64}}(undef, 0)
    Κ       = Vector{Vector{Float64}}(undef, 0)

    u = [Matrix(reshape(vec, 2, Int(length(vec) / 2))) for vec in sol.u]

    for ii in 1:s
        Area[ii] = Ω(u[ii])
        if Domain.domain_type == "2D"
            Cell_Count[ii] = size(u[ii], 2) / Domain.m
            Fnet, nV, den, stre, kap = PostCalcs2D(u[ii], p, AllCellMech[ii])
        else
            Cell_Count[ii] = (size(u[ii], 2) - 1) / Domain.m
            Fnet, nV, den, stre, kap = PostCalcs1D(u[ii], p, AllCellMech[ii])
        end
        push!(∑F, Fnet)
        push!(vₙ_out, nV)
        push!(DENSITY, vec(den))
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_t(Domain.btype, sol.t[1:s], u[1:s], ∑F, DENSITY, vₙ_out, Area, ψ, Κ, Cell_Count)
end


function postSimulation1D(sol, p, AllCellMech)
    Domain  = p[1]
    SimTime = p[3]

    c = size(sol.t, 1)
    s = min(c, size(AllCellMech, 1))
    Area       = Vector{Float64}(undef, s)
    Cell_Count = Vector{Float64}(undef, s)
    ∑F      = Vector{Vector{Float64}}(undef, 0)
    ψ       = Vector{Vector{Float64}}(undef, 0)
    DENSITY = Vector{Vector{Float64}}(undef, 0)
    vₙ_out  = Vector{Vector{Float64}}(undef, 0)
    Κ       = Vector{Vector{Float64}}(undef, 0)

    u = [Matrix(reshape(vec, 2, Int(length(vec) / 2))) for vec in sol.u]

    for ii in 1:s
        Cell_Count[ii] = (size(u[ii], 2) - 1) / Domain.m
        if Domain.domain_type == "2D"
            Fnet, nV, den, stre, kap = PostCalcs2D(u[ii], p, AllCellMech[ii])
        else
            Fnet, nV, den, stre, kap = PostCalcs1D(u[ii], p)
        end
        push!(∑F, Fnet)
        push!(vₙ_out, nV)
        push!(DENSITY, den)
        push!(ψ, stre)
        push!(Κ, kap)
    end
    return SimResults_t(Domain.btype, sol.t[1:s], u[1:s], ∑F, DENSITY, vₙ_out, Area, ψ, Κ, Cell_Count)
end


function calc_cell_orientation(embedded_cells)
    cell_orientation = Float64[]
    for cell in embedded_cells
        cell_right_edge = cell[:, end]
        cell_left_edge  = cell[:, 1]
        cell_from_origin = cell_right_edge - cell_left_edge
        cell_length = sqrt(sum((cell_right_edge - cell_left_edge).^2))
        θ_x = acos(clamp((cell_right_edge[1] - cell_left_edge[1]) / cell_length, -1.0, 1.0))

        if cell_from_origin[2] > 0
            ϕ_x = θ_x
        elseif cell_from_origin[2] < 0
            ϕ_x = 2π - θ_x
        else
            ϕ_x = cell_from_origin[1] > 0 ? 0.0 : π
        end
        isnan(ϕ_x) && error("Encountered NaN value for ϕ_x. Check the input data or calculations.")
        push!(cell_orientation, ϕ_x)
    end
    cell_orientation_out = round.(cell_orientation, digits=2)

    unique_cell_orientation = unique(cell_orientation_out)
    freq = [count(==(val), cell_orientation_out) for val in unique_cell_orientation]

    θ_care_range = collect(0:deg2rad(15):2π)
    θ_care_range .-= θ_care_range[2] / 2
    θ_care_range[1] = 2π + θ_care_range[1]

    freq_by_θ_care_range = Int[]
    for ii in 2:length(θ_care_range)
        v = 0
        for jj in eachindex(unique_cell_orientation)
            lo = ii == 2 ? θ_care_range[ii-1] - 2π : θ_care_range[ii-1]
            if lo <= unique_cell_orientation[jj] < θ_care_range[ii]
                v += freq[jj]
            end
        end
        push!(freq_by_θ_care_range, v)
    end

    return θ_care_range[2:end], freq_by_θ_care_range, cell_orientation_out
end


function calc_cell_orientation_at_t(embedded_cells, embed_cell_count, t_idx)
    cell_orientation = Float64[]
    n_cells = Int64(embed_cell_count[t_idx])
    for cell in embedded_cells[end-n_cells+1:end]
        cell_right_edge = cell[:, end]
        cell_left_edge  = cell[:, 1]
        cell_from_origin = cell_right_edge - cell_left_edge
        cell_length = sqrt(sum((cell_right_edge - cell_left_edge).^2))
        θ_x = acos(clamp((cell_right_edge[1] - cell_left_edge[1]) / cell_length, -1.0, 1.0))

        if cell_from_origin[2] > 0
            ϕ_x = θ_x
        elseif cell_from_origin[2] < 0
            ϕ_x = 2π - θ_x
        else
            ϕ_x = cell_from_origin[1] > 0 ? 0.0 : π
        end
        isnan(ϕ_x) && error("Encountered NaN value for ϕ_x. Check the input data or calculations.")
        push!(cell_orientation, ϕ_x)
    end
    cell_orientation_out = round.(cell_orientation, digits=2)

    unique_cell_orientation = unique(cell_orientation_out)
    freq = [count(==(val), cell_orientation_out) for val in unique_cell_orientation]

    θ_care_range = collect(0:deg2rad(15):2π)
    θ_care_range .-= θ_care_range[2] / 2
    θ_care_range[1] = 2π + θ_care_range[1]

    freq_by_θ_care_range = Int[]
    for ii in 2:length(θ_care_range)
        v = 0
        for jj in eachindex(unique_cell_orientation)
            lo = ii == 2 ? θ_care_range[ii-1] - 2π : θ_care_range[ii-1]
            if lo <= unique_cell_orientation[jj] < θ_care_range[ii]
                v += freq[jj]
            end
        end
        push!(freq_by_θ_care_range, v)
    end

    for r in freq_by_θ_care_range
        isnan(r) && error("Encountered NaN value for R. Check the input data or calculations.")
    end

    return θ_care_range[2:end], freq_by_θ_care_range, cell_orientation_out
end


function convert_coordinates_to_tuples(u)
    return [(u[1, i], u[2, i]) for i in axes(u, 2)]
end
