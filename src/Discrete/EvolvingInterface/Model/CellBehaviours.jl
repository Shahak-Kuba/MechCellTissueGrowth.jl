

function Set_Random_Seed(seednum=123)
    Random.seed!(seednum)
end

"""
    calc_spring_densities(uᵢ)

Calculate the spring densities for a given state vector `uᵢ`.

This function computes the inverse of the distance between each pair of adjacent elements in `uᵢ`. The calculation involves shifting `uᵢ` and then applying the δ function to each pair.

# Arguments
- `uᵢ`: A state vector representing the positions of particles or cells.

# Returns
A vector of spring densities, where each element is the inverse of the distance between adjacent elements in `uᵢ`.
"""
function calc_spring_densities(uᵢ)
    uᵢ₊₁ = circshift(uᵢ',1)
    return 1 ./ δ(uᵢ₊₁, uᵢ')
end

function calc_spring_lengths(uᵢ)
    uᵢ₊₁ = circshift(uᵢ',1)
    return δ(uᵢ₊₁, uᵢ')
end

#function calc_cell_lengths(u,m)
#    array = calc_spring_lengths(u)
#    return [sum(array[i:i+m-1]) for i in 1:m:length(array)-m+1]
#end

function calc_cell_lengths(u,m)
    spring_lengths = .√(sum((circshift(u,(0,-1))- u).^2,dims=1))
    return [sum(spring_lengths[i:i+m-1]) for i in 1:m:length(spring_lengths)-m+1]
end
"""
    calc_cell_densities(u, m)

Calculate the cell densities over a window of size `m` for a given state vector `u`.

This function first computes the spring densities using `calc_spring_densities` and then calculates the sum of these densities over a sliding window of size `m`. 

# Arguments
- `u`: A state vector representing the positions of particles or cells.
- `m`: The size of the window over which to calculate the densities.

# Returns
A vector of cell densities, each computed over a window of size `m`.
"""
function calc_cell_densities(u,m)
    # array = calc_spring_densities(u)
    return 1 ./ calc_cell_lengths(u,m)
end


"""
    P(event_f, ρ, α, fncs_type)

Calculate the probability of proliferation occurring, given the density `ρ`, the parameter `α`, and the function type `fncs_type`.

If the proliferation event is considered to occur (`event_f` is `true`), the probability is calculated based on the specified function type:
- For `"Constant"`: the probability is a constant value `α`.
- For `"Linear"`: the probability is scaled linearly by density `ρ`.
- For `"Nonlinear"`: the probability is inversely proportional to density `ρ`.

If the event is not considered to occur, a vector of zeros is returned.

### Arguments
- `event_f::Bool`: A boolean indicating whether the proliferation event is considered to occur.
- `ρ::AbstractVector`: A vector representing densities.
- `α::Float64`: Parameter for scaling the proliferation probability.
- `fncs_type::String`: Type of function used to calculate the probability. Options include:
  - `"Constant"`
  - `"Linear"`
  - `"Nonlinear"`

### Returns
A vector representing the calculated probabilities for proliferation.
"""
function P(event_f, ρ, α, fncs_type)
    if event_f
        if fncs_type == "Constant"
                return α.*ones(size(ρ))

        elseif fncs_type == "Constant2"
                lmin = 10.1
                P_rate = zeros(size(ρ))
                for cell in eachindex(ρ)
                    if 1 ./ ρ[cell] >= lmin
                        P_rate[cell] = α
                    end
                end
                return P_rate

        elseif fncs_type == "Length"
                return α.* (1 ./ ρ)
        end
    else
        return zeros(size(ρ))
    end
end


"""
    A(event_f, ρ, β, fncs_type)

Calculate the probability of apoptosis (cell death) occurring, given the density `ρ`, the parameter `β`, and the function type `fncs_type`.

If the apoptosis event is considered to occur (`event_f` is `true`), the probability is calculated based on the specified function type:
- For `"Constant"`: the probability is a constant value `β`.
- For `"Linear"`: the probability is scaled linearly by density `ρ`.
- For `"Nonlinear"`: the probability is inversely proportional to density `ρ`.

If the event is not considered to occur, a vector of zeros is returned.

### Arguments
- `event_f::Bool`: A boolean indicating whether the apoptosis event is considered to occur.
- `ρ::AbstractVector`: A vector representing densities.
- `β::Float64`: Parameter for scaling the apoptosis probability.
- `fncs_type::String`: Type of function used to calculate the probability. Options include:
  - `"Constant"`
  - `"Linear"`
  - `"Nonlinear"`

### Returns
A vector representing the calculated probabilities for apoptosis.
"""
function A(event_f, ρ, β, fncs_type)
    if event_f
       if fncs_type == "Constant"
            return β.*ones(size(ρ))
        elseif fncs_type == "Constant2"
            lmax = 10.1
            A_rate = zeros(size(ρ))
            for cell in eachindex(ρ)
                if 1 ./ ρ[cell] <= lmax
                    A_rate[cell] = β
                end
            end
            return A_rate
        elseif fncs_type == "Length"
            ld = 10.1
            A_rate = zeros(size(ρ))
            for ii in eachindex(ρ)
                if 1 ./ ρ[ii] <= ld
                    A_rate[ii] = β*(ld - (1 ./ ρ[ii]))
                else
                    A_rate[ii] = 0
                end
            end
            return A_rate
        end
    else
        return zeros(size(ρ))
    end
end


"""
    E(event_f, ρ, γ, fncs_type)

Calculate the probability of embedding occurring, given the density `ρ`, the parameter `γ`, and the function type `fncs_type`.

If the embedding event is considered to occur (`event_f` is `true`), the probability is calculated based on the specified function type:
- For `"Constant"`: the probability is a constant value `γ`.
- For `"Linear"`: the probability is scaled linearly by density `ρ`.
- For `"Nonlinear"`: the probability is inversely proportional to density `ρ`.

If the event is not considered to occur, a vector of zeros is returned.

### Arguments
- `event_f::Bool`: A boolean indicating whether the embedding event is considered to occur.
- `ρ::AbstractVector`: A vector representing cell densities.
- `γ::Float64`: Parameter for scaling the embedding probability.
- `fncs_type::String`: Type of function used to calculate the probability as functions of cell length l = 1/q. Options include:
  - `"Constant"`
  - `"Linear"`
  - `"Nonlinear"`

### Returns
A vector representing the calculated probabilities for embedding.
"""
function E(event_f, ρ, γ, fncs_type)
    if event_f
        if fncs_type == "Constant"
                return γ.*ones(size(ρ))
        elseif fncs_type == "Constant2"
                lmin = 0.0
                lmax = 16.1
                E_rate = zeros(size(ρ))
                for cell in eachindex(ρ)
                    if 1 ./ ρ[cell] <= lmax && 1 ./ ρ[cell] >= lmin
                        E_rate[cell] = γ
                    end
                end
                return E_rate
        elseif fncs_type == "Length"
            le_max = 20.0
            le_min = 10.0
            E_rate = zeros(size(ρ))
            for ii in eachindex(ρ)
                if 1 ./ ρ[ii] <= le_max 
                    #E_rate[ii] = γ*((le_max - le_min) - (1 ./ ρ[ii]))
                    E_rate[ii] = γ
                else
                    E_rate[ii] = 0.0
                end
            end
            return E_rate

        end
    else
        return zeros(size(ρ))
    end
end

function PE(event_f, ρ, γ, fncs_type)
    if event_f
        if fncs_type == "Constant"
            return γ.*ones(size(ρ))
        elseif fncs_type == "Constant2"
            lmin = 0.0
            lmax = 15.1
            PE_rate = zeros(size(ρ))
            for cell in eachindex(ρ)
                if 1 ./ ρ[cell] <= lmax && 1 ./ ρ[cell] >= lmin
                    PE_rate[cell] = γ
                end
            end
            return PE_rate
        end
    else
        return zeros(size(ρ))
    end
end
"""
    cell_probs(uᵢ, curr_t, m, SimTime, prolif, death, embed)

Calculate the probabilities for cell-related events (proliferation, death, embedding) based on the current state vector `uᵢ`.

This function uses `calc_cell_densities` to calculate cell densities and then computes the probabilities for proliferation, death, and embedding events. The probabilities are scaled by the time step size `δt`.

### Arguments
- `uᵢ::AbstractVector`: A state vector representing the positions of particles or cells.
- `curr_t::Float64`: The current simulation time.
- `m::Int64`: The size of the window over which to calculate the densities.
- `SimTime::SimTime_t`: An object containing simulation time parameters including event triggering and time step size.
- `prolif::CellEvent_t`: An object representing proliferation event parameters (flag, rate, and function type).
- `death::CellEvent_t`: An object representing death event parameters (flag, rate, and function type).
- `embed::CellEvent_t`: An object representing embedding event parameters (flag, rate, and function type).

### Returns
A tuple containing three vectors representing the calculated probabilities for:
1. Proliferation
2. Death
3. Embedding

If the respective event is not triggered based on the `event_trigger` setting, the function will return vectors of zeros for that event.
"""
function cell_probs(uᵢ,curr_t,m,SimTime,prolif,death,embed,prolifembed)
    ρ = calc_cell_densities(uᵢ,m)
    if SimTime.event_trigger == "Constant"
        return (P(prolif.flag,ρ,prolif.rate,prolif.event_func).*SimTime.event_δt, 
        A(death.flag,ρ,death.rate,death.event_func).*SimTime.event_δt, 
        E(embed.flag,ρ,embed.rate,embed.event_func).*SimTime.event_δt,
        PE(prolifembed.flag,ρ,prolifembed.rate,prolifembed.event_func).*SimTime.event_δt)
    elseif SimTime.event_trigger == "Periodic"
        if curr_t%SimTime.periodic_δt <= SimTime.event_length
            #println("Event executed at time: $curr_t")
            return (P(prolif.flag,ρ,prolif.rate,prolif.event_func).*SimTime.event_δt, 
                A(death.flag,ρ,death.rate,death.event_func).*SimTime.event_δt, 
                E(embed.flag,ρ,embed.rate,embed.event_func).*SimTime.event_δt,
                PE(prolifembed.flag,ρ,prolifembed.rate,prolifembed.event_func).*SimTime.event_δt)
        else
            p = zeros(length(P(prolif.flag,ρ,prolif.rate,prolif.event_func)))
            return (p, p, p, p)
        end
    end
end


"""
    find_cell_index(arr::Vector{Float64}, threshold::Float64)

Find the index in `arr` where the cumulative sum first exceeds or equals `threshold`.

This function is typically used in stochastic processes to determine an outcome based on a probability distribution.

# Arguments
- `arr`: A vector of probabilities.
- `threshold`: A threshold value used to find the corresponding index in `arr`.

# Returns
The index in `arr` where the cumulative sum first exceeds or equals `threshold`. If the threshold is not met, the length of `arr` is returned.

"""
function find_cell_index(arr::Vector{Float64}, threshold::Float64)
    cum_sum = cumsum(arr)
    index = findfirst(cum_sum .>= threshold)
    if index === nothing || index > length(arr)
        return index = length(arr)-1  # If the cumulative sum never reaches the threshold
    end
    return index
end


"""
    store_embed_cell_pos(pos)

Stores the position of embedded cells into embedded_cells array.

This function inserts the position vector `pos` of a boundary of an embedded cell into the embedded_cells array.

# Arguments
- `pos`: position vector of the cell that is being embedded into the tissue should contain `m+1` values given `m` springs in the simulation

# Returns
`nothing`. The function modifies `embedded_cells` in place.
"""
function store_embed_cell_pos(pos)
    global embedded_cells
    insert!(embedded_cells,size(embedded_cells,2),pos)
    return nothing
end

function store_embedded_cell(u, idx, m, t)
    global cell_embedment_times
    push!(cell_embedment_times, t)
    for i = idx:idx+m
        if i > size(u,2)
            store_embed_cell_pos(u[:,1].data)
        else
            store_embed_cell_pos(u[:,i].data)
        end
    end

    return nothing
end

function store_embedded_cell_count(u, t, integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = integrator.p
    global embedded_cells
    #global embedded_cell_count
    cell_count = size(hcat(embedded_cells...),2)/(Domain.m+1)
    #push!(embedded_cell_count,cell_count)
    return cell_count
end

function store_embed_rates(u,t,integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed= integrator.p
    u = integrator.u
    curr_t = integrator.t
    (p, a, e) = cell_probs(u, curr_t, Domain.m, SimTime, Prolif, Death, Embed, ProlifEmbed)
    return sum(e)/(length(e) * SimTime.event_δt)
end

function store_CellMech(u,t,integrator)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = integrator.p
    data = []
    push!(data,CellMech.kₛ)
    push!(data,CellMech.a)
    push!(data,CellMech.kf)
    push!(data,CellMech.growth_dir)
    return data
end

function convert_matrix(matrix, M)
    # Check if N is divisible by M
    N = size(matrix, 2)
    if N % M == 0 "N must be divisible by M"

        # Reshape the matrix
        reshaped_matrix = reshape(matrix, 2, M, :)

        # Split the reshaped matrix along the third dimension to get a vector of 2x3 matrices
        vector_of_matrices = [reshaped_matrix[:, :, i] for i in axes(reshaped_matrix,3)]

        return vector_of_matrices

    else
        return []
    end
end

function calc_cc_mechanism_rates(α,β,E_total)
    if α < 0 || α > 1
        error("α must be between 0 and 1")
    end
    if β < -α*E_total
        error("β must be greater than or equal to -α*Eᵗ")
    end
    P = β + α*E_total
    E = α*E_total
    PE = (1-α)*E_total
    return P, E, PE
end

