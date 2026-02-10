"""
    PostCalcs1D(u, p)

Perform post-calculation for 1D simulation data.

This function computes various physical quantities like force, density, velocity, stress, and curvature for a given state `u` and parameters `p`.

# Arguments
- `u`: A state vector representing the positions of particles or cells.
- `p`: A tuple of parameters used in the calculations.

# Returns
A tuple containing the sum of forces, normal velocity, density, stress, and curvature for each element in the state vector.
"""
function PostCalcs1D(u, p)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, aₘ, kₘ, ηₘ, kfₘ = p

    if Domain.btype == "InvertedBellCurve"
        dom = 1500; # For Bell curve
    else
        dom = 2*pi; # FOR Cosine SineWave
    end

    ∑F = zeros(size(u, 1))
    density = zeros(size(u, 1))
    vₙ = zeros(size(u, 1))
    ψ = zeros(size(u, 1))
    Κ = zeros(size(u, 1))

    uᵢ₊₁ = circshift(u,1)
    uᵢ₋₁ = circshift(u,-1)

    uᵢ₋₁[end,:] .= uᵢ₋₁[end,:] + [dom,0]
    uᵢ₊₁[1,:] .= uᵢ₊₁[1,:] - [dom,0]

    ∑F = diag((Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₘ,aₘ,CellMech.restoring_force))* transpose(τ(uᵢ₊₁,uᵢ₋₁)))
    #diag(((Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force)) + (Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force)) )* transpose(τ(uᵢ₊₁,uᵢ₋₁)))
    density = (ρ(uᵢ₊₁, u).+ρ(u, uᵢ₋₁))./(2*Domain.m)
    ψ = ∑F / (kₘ*aₘ)
    Κ = κ(uᵢ₋₁,u,uᵢ₊₁)
    vₙx = Vₙ(uᵢ₋₁,u,uᵢ₊₁,kfₘ,SimTime.δt,"2D")[:,1]
    vₙy = Vₙ(uᵢ₋₁,u,uᵢ₊₁,kfₘ,SimTime.δt,"2D")[:,2]
    vₙ = .√(vₙx.^2 + vₙy.^2)
    

    return ∑F, vₙ, density, ψ, Κ
end


###################################################################################################

"""
    PostCalcs2D(u, p)

Perform post-calculation for 2D simulation data.

This function is similar to `PostCalcs1D`, but it is tailored for 2D simulation data. It calculates force, density, velocity, stress, and curvature in a 2D context.

# Arguments
- `u`: A 2D state vector representing the positions of particles or cells.
- `p`: A tuple of parameters used in the calculations.

# Returns
A tuple containing the sum of forces, normal velocity, density, stress, and curvature for each element in the state vector.
"""
function PostCalcs2D(u, p, CellMech_at_t)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = p

    #u = reshape(u, Int(length(u)/2), 2)
    kₘ = CellMech_at_t[1]
    aₘ = CellMech_at_t[2]
    kfₘ = CellMech_at_t[3]
    growth_dir = CellMech_at_t[4]
    η = CellMech.η

    ∑F = zeros(size(u, 2))
    density = zeros(size(u, 2))
    vₙ = zeros(size(u, 2))
    ψ = zeros(size(u, 2))
    Κ = zeros(size(u, 2))

    uᵢ₊₁ = circshift(u,(0,-1))
    uᵢ₋₁ = circshift(u,(0,1))

    l,ρ,dv,τ,n = calc_l_ρ_dv_τ_n(uᵢ₊₁, u, uᵢ₋₁, growth_dir)

    ∑F .= diag(( Fₛ⁺(l,kₘ,aₘ,CellMech.restoring_force).* dv + Fₛ⁻(circshift(l,(0,1)),kₘ,aₘ,CellMech.restoring_force).*circshift(dv,(0,1)) )' * τ)
    #diag(((Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force)) + (Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force)) )* transpose(τ(uᵢ₊₁,uᵢ₋₁)))
    density .= ((circshift(ρ,(0, -1)).+ρ)/(2*Domain.m))'
    ψ .= ∑F ./ (kₘ.*aₘ)
    Κ .= κ(uᵢ₋₁,u,uᵢ₊₁)
    vₙx = Vₙ(uᵢ₋₁,u,uᵢ₊₁,ρ,n,kfₘ,SimTime.δt)[1,:]
    vₙy = Vₙ(uᵢ₋₁,u,uᵢ₊₁,ρ,n,kfₘ,SimTime.δt)[2,:]
    vₙ .= .√(vₙx.^2 + vₙy.^2)

    return ∑F, vₙ, density, ψ, Κ
end

"""
    postSimulation2D(btype, sol, p)

Perform post simulation calculations for 2D simulation and return a comprehensive data structure with all relevant data.

This function processes the solution from a 2D simulation, similarly to `postSimulation1D`, but adapted for 2D data.

# Arguments
- `btype`: The type of boundary condition or simulation.
- `sol`: The solution object from the simulation.
- `p`: Parameters used in the post calculations.

# Returns
An instance of `SimResults_t` containing the calculated data.
"""
function postSimulation(sol, p, AllCellMech)

    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = p

    c = size(sol.t, 1)
    s = min(c, size(AllCellMech,1))
    Area = Vector{Float64}(undef, s)
    Cell_Count = Vector{Float64}(undef, s)
    ∑F = Vector{Vector{Float64}}(undef, 0)
    ψ = Vector{Vector{Float64}}(undef, 0)
    DENSITY = Vector{Vector{Float64}}(undef, 0)
    vₙ = Vector{Vector{Float64}}(undef, 0)
    Κ = Vector{Vector{Float64}}(undef, 0)

    u = [Matrix((reshape(vec, 2, Int(length(vec)/2)))) for vec in sol.u]

    # adding periodic boundary node in 1D case
    if Domain.domain_type == "1D"
        if Domain.btype == "InvertedBellCurve"
            dom = 1500
        else
            dom = 2π
        end
        u = [[vec; (vec[1,:] + [dom,0])'] for vec in u]
    end


    for ii in 1:s
        Area[ii] = Ω(u[ii]) # area calculation
        Cell_Count[ii] = size(u[ii],2)/Domain.m
        if Domain.domain_type == "2D"
            Fnet, nV, den, stre, kap = PostCalcs2D(u[ii], p, AllCellMech[ii])
        else
            Fnet, nV, den, stre, kap = PostCalcs1D(u[ii], p)
        end

        push!(∑F, Fnet)
        push!(vₙ, nV)
        push!(DENSITY, den)
        push!(ψ, stre)
        push!(Κ, kap)
    end

    return SimResults_t(Domain.btype, sol.t[1:s], u[1:s], ∑F, DENSITY, vₙ, Area, ψ, Κ, Cell_Count)
end


function calc_cell_orientation(embedded_cells)
    cell_orientation = Float64[]
    for cell in embedded_cells
        cell_right_edge = cell[:,end]
        cell_left_edge = cell[:,1]
        cell_from_origin = cell_right_edge - cell_left_edge
        cell_length = .√(sum((cell_right_edge - cell_left_edge).^2, dims=1))[1] 
        θ_x = acos((cell_right_edge[1] - cell_left_edge[1]) / cell_length)
        # checking what quadrant the cell orientation is in
        ϕ_x = 0
        if cell_from_origin[2] > 0
            ϕ_x = θ_x
        elseif cell_from_origin[2] < 0
            ϕ_x = 2π - θ_x
        else
            if cell_from_origin[1] > 0
                ϕ_x = 0
            else
                ϕ_x = π
            end
        end
        if isnan(ϕ_x)
            error("Encountered NaN value for ϕ_x. Check the input data or calculations.")
        end
        push!(cell_orientation, ϕ_x)
    end
    cell_orientation_out = round.(cell_orientation, digits=2)

    unique_cell_orientation = unique(round.(cell_orientation, digits=2))
    freq = Int64[]
    for val in unique_cell_orientation
        push!(freq, count(i -> (i == val), round.(cell_orientation, digits=2)))
    end

    θ_care_range = 0:deg2rad(15):2π
    θ_care_range = hcat([θ_care_range .- θ_care_range[2]/2]...)
    θ_care_range[1] = 2π + θ_care_range[1]

    freq_by_θ_care_range = []
    for ii in 2:size(θ_care_range, 1)
        v = 0
        for jj in eachindex(unique_cell_orientation)
            if ii == 2
                if (θ_care_range[ii-1] - 2π <= unique_cell_orientation[jj]) && (unique_cell_orientation[jj] < θ_care_range[ii])
                    v += freq[jj]
                end
            else
                if (θ_care_range[ii-1] <= unique_cell_orientation[jj]) && (unique_cell_orientation[jj] < θ_care_range[ii])
                    v += freq[jj]
                end
            end
        end
        push!(freq_by_θ_care_range, v)
    end

    θ = vcat(θ_care_range)[2:end]
    R = vcat(freq_by_θ_care_range)

    return θ, R, cell_orientation_out
end

function calc_cell_orientation_at_t(embedded_cells, embed_cell_count, t_idx)
    cell_orientation = Float64[]
    #println(size(embedded_cells))
    for cell in embedded_cells[end - Int64((embed_cell_count[t_idx])[1]) + 1:end]
        cell_right_edge = cell[:,end]
        cell_left_edge = cell[:,1]
        cell_from_origin = cell_right_edge - cell_left_edge
        cell_length = .√(sum((cell_right_edge - cell_left_edge).^2, dims=1))[1] 
        θ_x = acos((cell_right_edge[1] - cell_left_edge[1]) / cell_length)
        # checking what quadrant the cell orientation is in
        ϕ_x = 0
        if cell_from_origin[2] > 0
            ϕ_x = θ_x
        elseif cell_from_origin[2] < 0
            ϕ_x = 2π - θ_x
        else
            if cell_from_origin[1] > 0
                ϕ_x = 0
            else
                ϕ_x = π
            end
        end
        if isnan(ϕ_x)
            error("Encountered NaN value for ϕ_x. Check the input data or calculations.")
        end
        push!(cell_orientation, ϕ_x)
    end
    cell_orientation_out = round.(cell_orientation, digits=2)

    unique_cell_orientation = unique(round.(cell_orientation, digits=2))
    freq = Int64[]
    for val in unique_cell_orientation
        push!(freq, count(i -> (i == val), round.(cell_orientation, digits=2)))
    end

    θ_care_range = 0:deg2rad(15):2π
    θ_care_range = hcat([θ_care_range .- θ_care_range[2]/2]...)
    θ_care_range[1] = 2π + θ_care_range[1]

    freq_by_θ_care_range = []
    for ii in 2:size(θ_care_range, 1)
        v = 0
        for jj in eachindex(unique_cell_orientation)
            if ii == 2
                if (θ_care_range[ii-1] - 2π <= unique_cell_orientation[jj]) && (unique_cell_orientation[jj] < θ_care_range[ii])
                    v += freq[jj]
                end
            else
                if (θ_care_range[ii-1] <= unique_cell_orientation[jj]) && (unique_cell_orientation[jj] < θ_care_range[ii])
                    v += freq[jj]
                end
            end
        end
        push!(freq_by_θ_care_range, v)
    end

    θ = vcat(θ_care_range)[2:end]
    R = vcat(freq_by_θ_care_range)
    for ii in eachindex(R)
        if isnan(R[ii])
            error("Encountered NaN value for R. Check the input data or calculations.")
        end
    end

    return θ, R, cell_orientation_out   
end

function convert_coordinates_to_tuples(u)
    coordinates = [(u[1, i], u[2, i]) for i in axes(u, 2)]
    return coordinates
end