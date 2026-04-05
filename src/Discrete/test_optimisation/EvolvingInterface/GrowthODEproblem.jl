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
function Growth_ODE!(du, u, p, t)
    Domain, CellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, cache, state = p
    N   = size(u, 2)
    η_inv = 1.0 / CellMech.η
    δt  = SimTime.δt

    # Resize cache if the number of nodes changed (after cell events)
    if cache.N != N
        resize_cache!(cache, N)
    end

    if Domain.domain_type == "2D"
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

            # Normal vector = 90° rotation of dv, sign from growth_dir
            gs = CellMech.growth_dir[j]          # Int8: -1 or +1
            cache.nx[j] =  gs * cache.dvy[j]
            cache.ny[j] = -gs * cache.dvx[j]
        end

        # ── Pass 2: Compute du = mechanical + growth velocity ────────
        @inbounds for j in 1:N
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

            # --- Growth (normal) velocity via line intersection ---
            # Left spring (j-1 → j)
            ρ_left  = cache.ρ[jm1]
            V_left  = CellMech.kf[j] * ρ_left
            nl_x    = cache.nx[jm1]
            nl_y    = cache.ny[jm1]

            # Right spring (j → j+1)
            ρ_right = cache.ρ[j]
            V_right = CellMech.kf[j] * ρ_right
            nr_x    = cache.nx[j]
            nr_y    = cache.ny[j]

            # Offset points for line intersection
            vl_nx = V_left  * nl_x * δt
            vl_ny = V_left  * nl_y * δt
            vr_nx = V_right * nr_x * δt
            vr_ny = V_right * nr_y * δt

            rm1x = u[1, j]   + vl_nx
            rm1y = u[2, j]   + vl_ny
            rlx  = u[1, jm1] + vl_nx
            rly  = u[2, jm1] + vl_ny

            rm2x = u[1, j]   + vr_nx
            rm2y = u[2, j]   + vr_ny
            rrx  = u[1, jp1] + vr_nx
            rry  = u[2, jp1] + vr_ny

            ix, iy = _line_intersect(rm1x, rm1y, rlx, rly,
                                     rm2x, rm2y, rrx, rry)

            growth_vx = (ix - u[1, j]) / δt
            growth_vy = (iy - u[2, j]) / δt

            # --- Total ---
            du[1, j] = mech_vx + growth_vx
            du[2, j] = mech_vy + growth_vy
        end
    end
    nothing
end
