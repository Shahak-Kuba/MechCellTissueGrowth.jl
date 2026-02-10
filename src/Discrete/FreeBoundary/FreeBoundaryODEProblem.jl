function FB_ODE!(du, u, p , t)
    TissueMech, SimTime, y = p
    uᵢ₊₁ = circshift(u,(0,-1))
    uᵢ₋₁ = circshift(u,(0,1))

    l,ρ,dv,τ,n = calc_l_ρ_dv_τ_n(uᵢ₊₁, u, uᵢ₋₁, "inward")

    du[:,1] = [0. ; 0.]
    du[:,2:end-1] = (1 ./ TissueMech.η[1]) .* ( -F(l[1:end-2],TissueMech.k[1:end-2],TissueMech.a[1:end-2],TissueMech.restoring_force).* dv[:,1:end-2]       +     F(l[2:end-1],TissueMech.k[2:end-1],TissueMech.a[2:end-1],TissueMech.restoring_force).* dv[:,2:end-1] )
    du[:,end] = (1 ./ TissueMech.η[1]) .* ( -F(l[end-1],TissueMech.k[end-1],TissueMech.a[end-1],TissueMech.restoring_force).* dv[:,end-1] )
    
end