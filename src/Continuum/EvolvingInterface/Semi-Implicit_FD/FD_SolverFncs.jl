function aₘ(Rᵢ₊₁, Rᵢ₋₁)
    a₊ = zeros(size(Rᵢ₊₁))
    a₋ = zeros(size(Rᵢ₊₁))
    for i in eachindex(Rᵢ₊₁)
        if Rᵢ₊₁[i] - Rᵢ₋₁[i] > 0
            a₊[i] = 1.0;
            a₋[i] = 0.0;
        else
            a₊[i] = 0.0;
            a₋[i] = 1.0;
        end
    end
    return a₊,  a₋
end

# same function created in FVM_SolverFncs.jl

function FD_InitialBoundary(type,r₀,θ,M)
    r = zeros(1,M);

    if type == "circle"
        r .= ones(size(r)).*r₀
    elseif type == "square"
        R = 50;
        for i in 1:M
            if θ[i] < π/2
                r[i] = min((R/cos(θ[i])), (R/sin(θ[i])))
            else
                r[i] = min((R/cos((θ[i]%(π/2)))), (R/sin((θ[i]%(π/2)))))
            end
        end
    elseif type =="hex"
        R = √((2/3√3)*π*(r₀^2));
        for i in 1:M
            if θ[i] >= 0 && θ[i] < π/3
                r[i] = R/(cos(θ[i])+1/sqrt(3)*sin(θ[i]))
            elseif θ[i] >= π/3 && θ[i] < π/2
                r[i] = R*sqrt(3)/(2*sin(θ[i]))
            elseif θ[i] >= π/2 && θ[i] < 2*π/3
                r[i] = R*sqrt(3)/(2*sin(θ[i]))
            elseif θ[i] >= 2*π/3 && θ[i] < π
                r[i] = R/(-cos(θ[i])+1/sqrt(3)*sin(θ[i]))
            elseif θ[i] >= π && θ[i] < 4*π/3
                r[i] = -R/(cos(θ[i])+1/sqrt(3)*sin(θ[i]))
            elseif θ[i] >= 4*π/3 && θ[i] < 3*π/2
                r[i] = -R*sqrt(3)/(2*sin(θ[i])) 
            elseif θ[i] >= 3*π/2 && θ[i] < 5*π/3
                r[i] = -R*sqrt(3)/(2*sin(θ[i]))
            else
                r[i] = R/(cos(θ[i])-1/sqrt(3)*sin(θ[i]))
            end
        end
                       
    end
    return r
end
