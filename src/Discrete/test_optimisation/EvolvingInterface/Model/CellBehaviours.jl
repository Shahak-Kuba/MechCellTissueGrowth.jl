

function Set_Random_Seed(seednum=123)
    Random.seed!(seednum)
end

"""
    calc_spring_densities(uᵢ)

Calculate the spring densities (1/length) for a given state matrix `uᵢ`.
"""
function calc_spring_densities(uᵢ)
    uᵢ₊₁ = circshift(uᵢ', 1)
    return 1.0 ./ δ(uᵢ₊₁, uᵢ')
end

function calc_spring_lengths(uᵢ)
    uᵢ₊₁ = circshift(uᵢ', 1)
    return δ(uᵢ₊₁, uᵢ')
end

function calc_cell_lengths(u, m)
    spring_lengths = .√(sum((circshift(u, (0, -1)) .- u).^2, dims=1))
    return [sum(spring_lengths[i:i+m-1]) for i in 1:m:length(spring_lengths)-m+1]
end

"""
    calc_cell_densities(u, m)

Calculate cell densities over a window of size `m`.
"""
function calc_cell_densities(u, m)
    return 1.0 ./ calc_cell_lengths(u, m)
end


"""
    P(event_f, ρ, α, fncs_type)

Calculate proliferation probability given density `ρ`, rate `α`, and function type.
"""
function P(event_f, ρ, α, fncs_type::EventFunc)
    !event_f && return zeros(length(ρ))

    if fncs_type == EF_Constant
        return fill(α, length(ρ))
    elseif fncs_type == EF_Constant2
        lmin = 10.1
        P_rate = zeros(length(ρ))
        @inbounds for cell in eachindex(ρ)
            if 1.0 / ρ[cell] >= lmin
                P_rate[cell] = α
            end
        end
        return P_rate
    elseif fncs_type == EF_Length
        return α .* (1.0 ./ ρ)
    end
    return zeros(length(ρ))
end


"""
    A(event_f, ρ, β, fncs_type)

Calculate apoptosis (cell death) probability given density `ρ`, rate `β`, and function type.
"""
function A(event_f, ρ, β, fncs_type::EventFunc)
    !event_f && return zeros(length(ρ))

    if fncs_type == EF_Constant
        return fill(β, length(ρ))
    elseif fncs_type == EF_Constant2
        lmax = 10.1
        A_rate = zeros(length(ρ))
        @inbounds for cell in eachindex(ρ)
            if 1.0 / ρ[cell] <= lmax
                A_rate[cell] = β
            end
        end
        return A_rate
    elseif fncs_type == EF_Length
        ld = 10.1
        A_rate = zeros(length(ρ))
        @inbounds for ii in eachindex(ρ)
            if 1.0 / ρ[ii] <= ld
                A_rate[ii] = β * (ld - 1.0 / ρ[ii])
            end
        end
        return A_rate
    end
    return zeros(length(ρ))
end


"""
    E(event_f, ρ, γ, fncs_type)

Calculate embedding probability given density `ρ`, rate `γ`, and function type.
"""
function E(event_f, ρ, γ, fncs_type::EventFunc)
    !event_f && return zeros(length(ρ))

    if fncs_type == EF_Constant
        return fill(γ, length(ρ))
    elseif fncs_type == EF_Constant2
        lmin = 0.0
        lmax = 16.1
        E_rate = zeros(length(ρ))
        @inbounds for cell in eachindex(ρ)
            cell_len = 1.0 / ρ[cell]
            if cell_len <= lmax && cell_len >= lmin
                E_rate[cell] = γ
            end
        end
        return E_rate
    elseif fncs_type == EF_Length
        le_max = 20.0
        E_rate = zeros(length(ρ))
        @inbounds for ii in eachindex(ρ)
            if 1.0 / ρ[ii] <= le_max
                E_rate[ii] = γ
            end
        end
        return E_rate
    end
    return zeros(length(ρ))
end

function PE(event_f, ρ, γ, fncs_type::EventFunc)
    !event_f && return zeros(length(ρ))

    if fncs_type == EF_Constant
        return fill(γ, length(ρ))
    elseif fncs_type == EF_Constant2
        lmin = 0.0
        lmax = 15.1
        PE_rate = zeros(length(ρ))
        @inbounds for cell in eachindex(ρ)
            cell_len = 1.0 / ρ[cell]
            if cell_len <= lmax && cell_len >= lmin
                PE_rate[cell] = γ
            end
        end
        return PE_rate
    end
    return zeros(length(ρ))
end

"""
    cell_probs(uᵢ, curr_t, m, SimTime, prolif, death, embed, prolifembed)

Calculate probabilities for cell events (proliferation, death, embedding)
based on the current state vector `uᵢ`.
"""
function cell_probs(uᵢ, curr_t, m, SimTime, prolif, death, embed, prolifembed)
    ρ = calc_cell_densities(uᵢ, m)
    δt = SimTime.event_δt

    if SimTime.event_trigger == ET_Constant
        return (P(prolif.flag, ρ, prolif.rate, prolif.event_func) .* δt,
                A(death.flag, ρ, death.rate, death.event_func) .* δt,
                E(embed.flag, ρ, embed.rate, embed.event_func) .* δt,
                PE(prolifembed.flag, ρ, prolifembed.rate, prolifembed.event_func) .* δt)
    elseif SimTime.event_trigger == ET_Periodic
        if curr_t % SimTime.periodic_δt <= SimTime.event_length
            return (P(prolif.flag, ρ, prolif.rate, prolif.event_func) .* δt,
                    A(death.flag, ρ, death.rate, death.event_func) .* δt,
                    E(embed.flag, ρ, embed.rate, embed.event_func) .* δt,
                    PE(prolifembed.flag, ρ, prolifembed.rate, prolifembed.event_func) .* δt)
        else
            zp = zeros(length(ρ))
            return (zp, zp, zp, zp)
        end
    end
end


"""
    store_embed_cell_pos(state, pos)

Stores the position of an embedded cell node into the simulation state.
"""
function store_embed_cell_pos(state::SimulationState, pos)
    push!(state.embedded_cells, pos)
    return nothing
end

function store_embedded_cell(state::SimulationState, u, idx, m, t)
    push!(state.cell_embedment_times, t)
    for i = idx:idx+m
        if i > size(u, 2)
            store_embed_cell_pos(state, Vector(u[:, 1]))
        else
            store_embed_cell_pos(state, Vector(u[:, i]))
        end
    end
    return nothing
end

function store_embedded_cell_count(u, t, integrator)
    Domain = integrator.p[1]
    state  = integrator.p[9]
    if isempty(state.embedded_cells)
        return 0.0
    end
    cell_count = length(state.embedded_cells) / (Domain.m + 1)
    return cell_count
end

function store_embed_rates(u, t, integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = integrator.p[1:7]
    u_curr    = integrator.u
    curr_t    = integrator.t
    (p_prob, a_prob, e_prob, pe_prob) = cell_probs(u_curr, curr_t, Domain.m, SimTime, Prolif, Death, Embed, ProlifEmbed)
    return sum(e_prob) / (length(e_prob) * SimTime.event_δt)
end

function store_CellMech(u, t, integrator)
    CellMech = integrator.p[2]
    return (collect(CellMech.kₛ), collect(CellMech.a), collect(CellMech.kf), collect(CellMech.growth_dir))
end

function convert_matrix(matrix, M)
    N = size(matrix, 2)
    if N % M != 0
        return Matrix{Float64}[]
    end
    reshaped_matrix = reshape(matrix, 2, M, :)
    return [reshaped_matrix[:, :, i] for i in axes(reshaped_matrix, 3)]
end

function calc_cc_mechanism_rates(α, β, E_total)
    if α < 0 || α > 1
        error("α must be between 0 and 1")
    end
    if β < -α * E_total
        error("β must be greater than or equal to -α*Eᵗ")
    end
    P = β + α * E_total
    E = α * E_total
    PE = (1 - α) * E_total
    return P, E, PE
end
