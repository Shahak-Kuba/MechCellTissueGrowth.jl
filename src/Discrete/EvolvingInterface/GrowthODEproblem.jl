"""
    Growth_ODE!(du, u, p, t)

Define the ODE system for 2D mechanical relaxation with initial conditions. This function computes the derivatives `du` based on the current state `u` and parameters `p`.

# Arguments
- `du`: Array to store the derivatives of `u`.
- `u`: Current position array of all spring boundary nodes.
- `p`: Parameters tuple `(N, kₛ, η, kf, l₀, δt, growth_dir)`.
- `t`: Current time.

# Description
Calculates the mechanical relaxation and normal velocity in a 1D system with periodic boundary conditions. The derivatives are based on spring forces and mechanical properties defined in `p`.
"""

"""
function Growth_ODE!(du,u,p,t) 
    m,kₛ,η,kf,l₀,δt,growth_dir,domain_type,btype,restoring_force,prolif,death,embed = p
    uᵢ₊₁ = circshift(u',1)
    uᵢ₋₁ = circshift(u',-1)

    if domain_type == "2D"
        du .= ((1/η) .* diag((Fₛ⁺(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force) + Fₛ⁻(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁) +
                       Vₙ(uᵢ₋₁,u',uᵢ₊₁,kf,δt,growth_dir))'
        # Below is unstable normal vlocity
        #du .= ((1/η) .* diag((Fₛ⁺(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force) + Fₛ⁻(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁) +
        #                Vₙ(uᵢ₋₁, u', uᵢ₊₁, kf,"inward"))'
    else
        if btype == "InvertedBellCurve"
            dom = 1500; # For Bell curve
        else
            dom = 2*pi; # FOR Cosine SineWave
        end
        uᵢ₋₁[end,:] = uᵢ₋₁[end,:]+[dom;0]
        uᵢ₊₁[1,:] = uᵢ₊₁[1,:]-[dom;0]
        du .= ((1/η) .* diag((Fₛ⁺(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force) + Fₛ⁻(u',uᵢ₊₁,uᵢ₋₁,kₛ,l₀,restoring_force)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁) +
                            Vₙ(uᵢ₋₁,u',uᵢ₊₁,kf,δt,growth_dir))'
    end
    nothing
end
"""

function Growth_ODE!(du,u,p,t) 
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed = p
    uᵢ₊₁ = circshift(u,(0,-1))
    uᵢ₋₁ = circshift(u,(0,1))

    l,ρ,dv,τ,n = calc_l_ρ_dv_τ_n(uᵢ₊₁, u, uᵢ₋₁, CellMech.growth_dir)

    if Domain.domain_type == "2D"
        #du .= (((1/ηₘ) .* diag(( KubaPhD.Fₛ⁺(l,kₘ,aₘ,dv,CellMech.restoring_force) + KubaPhD.Fₛ⁻(circshift(l,(0,1)),kₘ,aₘ,circshift(dv,(0,1)),CellMech.restoring_force) )' * τ ))' .* τ) +
        #               Vₙ(uᵢ₋₁,u,uᵢ₊₁,ρ,n,kfₘ,SimTime.δt)
        
        du .= (1/CellMech.η) .* (( diag(( Fₛ⁺(l,CellMech.kₛ,CellMech.a,CellMech.restoring_force).* dv + Fₛ⁻(circshift(l,(0,1)),CellMech.kₛ,CellMech.a,CellMech.restoring_force).*circshift(dv,(0,1)) )' * τ) )' .* τ) .+  Vₙ(uᵢ₋₁,u,uᵢ₊₁,ρ,n,CellMech.kf,SimTime.δt)
        #du .= Vₙ(uᵢ₋₁,u,uᵢ₊₁,ρ,n,kfₘ,SimTime.δt)

    else # consider 1D (NEED TO DO)
        #if Domain.btype == "InvertedBellCurve"
        #    dom = 1500; # For Bell curve
        #else
        #    dom = 2*pi; # FOR Cosine SineWave
        #end
        #uᵢ₋₁[end,:] = uᵢ₋₁[end,:]+[dom;0]
        #uᵢ₊₁[1,:] = uᵢ₊₁[1,:]-[dom;0]
        #du .= ((1/ηₘ) .* diag((Fₛ⁺(u,uᵢ₊₁,uᵢ₋₁,kₘ,aₘ,CellMech.restoring_force) + Fₛ⁻(u,uᵢ₊₁,uᵢ₋₁,kₘ,aₘ,CellMech.restoring_force)) * transpose(τ(uᵢ₊₁,uᵢ₋₁))).*τ(uᵢ₊₁,uᵢ₋₁) +
        #               Vₙ(uᵢ₋₁,u,uᵢ₊₁,kfₘ,SimTime.δt,CellMech.growth_dir))'
    end
    nothing
end



