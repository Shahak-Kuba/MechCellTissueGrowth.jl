"""
function FB_ODE!(du, u, p, t)
    # p is a 9-element tuple — we only need the first 7 fields here.
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, cache, state = p

    uᵢ₊₁ = circshift(u, (0, -1))
    uᵢ₋₁ = circshift(u, (0,  1))

    l, ρ_val, dv, τ_val, n_val = calc_l_ρ_dv_τ_n(uᵢ₊₁, u, uᵢ₋₁, "inward")

    du[:, 1] = [0.0; 0.0]
    du[:, 2:end-1] = (1.0 ./ CellMech.η) .* (
        -F(l[1:end-2]', CellMech.kₛ[1:end-2], CellMech.a[1:end-2], CellMech.restoring_force) .* dv[:, 1:end-2] .+
         F(l[2:end-1]', CellMech.kₛ[2:end-1], CellMech.a[2:end-1], CellMech.restoring_force) .* dv[:, 2:end-1])
    du[:, end] = (1.0 ./ CellMech.η) .* (
        -F(l[end-1], CellMech.kₛ[end-1], CellMech.a[end-1], CellMech.restoring_force) .* dv[:, end-1])
end
"""

"""
    Growth_ODE!(du, u, p, t)

Zero-allocation ODE right-hand side for 2D mechanical relaxation with tissue growth.

All intermediate quantities are computed in-place using scalar loops and
the pre-allocated `ODECache` stored in `p`.

Key optimisations vs. the original:
- No `circshift` — uses manual circular indexing.
- No `diag(A' * B)` — computes column-wise dot products inline.
- No temporary arrays — all buffers live in `cache`.
- No string comparison — uses `ForceType` enum for force dispatch.
"""
function FB_ODE!(du, u, p, t)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, cache, state = p
    N   = size(u, 2)
    η_inv = 1.0 / CellMech.η
    δt  = SimTime.δt

    # Resize cache if the number of nodes changed (after cell events)
    if cache.N != N
        resize_cache!(cache, N)
    end

        # ── Pass 1: Compute spring quantities ────────────────────────
        @inbounds for j in 1:N
            jp1 = j == N ? 1 : j + 1
            jm1 = j == 1 ? N : j - 1

            # Spring length  j → j+1
            dx = u[1, jp1] - u[1, j]
            dy = u[2, jp1] - u[2, j]
            lj = sqrt(dx * dx + dy * dy)
            cache.l[j]   = lj
            cache.ρ[j]   = 1.0 / lj
            cache.dvx[j] = dx / lj
            cache.dvy[j] = dy / lj

            # Tangent vector (j-1 → j+1)
            tx = u[1, jp1] - u[1, jm1]
            ty = u[2, jp1] - u[2, jm1]
            tl = sqrt(tx * tx + ty * ty)
            cache.τx[j] = tx / tl
            cache.τy[j] = ty / tl
        end

        # ── Pass 2: Compute du = mechanical────────
        @inbounds for j in 1:N-1
            jp1 = j == N ? 1 : j + 1
            jm1 = j == 1 ? N : j - 1

            # --- Mechanical (tangential) velocity ---
            # Force from spring j  (right neighbour)
            F_plus  = force_scalar(cache.l[j],   CellMech.kₛ[j],   CellMech.a[j],   CellMech.restoring_force)
            # Force from spring j-1 (left neighbour)
            F_minus = force_scalar(cache.l[jm1], CellMech.kₛ[jm1], CellMech.a[jm1], CellMech.restoring_force)

            # Net force vector
            fx = F_plus * cache.dvx[j] - F_minus * cache.dvx[jm1]
            fy = F_plus * cache.dvy[j] - F_minus * cache.dvy[jm1]

            # Project onto tangent
            F_tang = fx * cache.τx[j] + fy * cache.τy[j]

            mech_vx = η_inv * F_tang * cache.τx[j]
            mech_vy = η_inv * F_tang * cache.τy[j]

            # --- Total ---
            du[1, j] = mech_vx
            du[2, j] = mech_vy
        end

        # for final node
        # --- Mechanical (tangential) velocity ---
        
            # Force from spring j-1 (left neighbour)
            F_minus = force_scalar(cache.l[N-1], CellMech.kₛ[N-1], CellMech.a[N-1], CellMech.restoring_force)

            # Net force vector
            fx = - F_minus * cache.dvx[N-1]
            fy = - F_minus * cache.dvy[N-1]

            # Project onto tangent
            F_tang = fx * cache.τx[N-1] + fy * cache.τy[N-1]

            mech_vx = η_inv * F_tang * cache.τx[N-1]
            mech_vy = η_inv * F_tang * cache.τy[N-1]

            # --- Total ---
            du[1, N] = mech_vx
            du[2, N] = mech_vy

    nothing
end

