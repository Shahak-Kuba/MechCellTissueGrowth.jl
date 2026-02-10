# Analytic solution for number of cells
N_analytic(N₀,P,A,E,t) = N₀.*exp.((P.-(A.+E)).*t)

# Analytic solution for void area
Ω_analytic(Ω₀,N₀,kf,t) = Ω₀ .- kf.*N₀.*t # constant cell population
Ω_analytic(Ω₀,N₀,A,kf,Ot,t) = Ω₀ .- (N₀./(A.+kf.*Ot)).*(1 .- exp.(-(A+kf*Ot).*t)) # changing cell population

# function to calculate kf value from experiment

Rr(KF, Ω₀, t)= sqrt.((Ω₀.-KF.*t)/π)

V(KF, Ω₀, t) = (-KF ./ (2 .*(π^(0.5))) .* ((1)./((Ω₀ .- KF.*t).^(0.5))))

V(KF, Ω₀, Tb, t) = (KF ./ 2π) .* ((Ω₀ - KF.*(Tb .- t))./(π)).^(0.5)

Ncells(KF, Ω₀, ρ₀, t) = 2 .* π .* ρ₀ .* ((Ω₀ .- KF.*t) ./ π).^(0.5)


function q_analytic(kf, q₀, L, Ω₀, P, A, E, Tb, t)
    q_t = zeros(length(t))
    for i in eachindex(t)
        if 0 <= 1 - (kf*L*q₀)/(Ω₀*(P - (A + E))) * (exp((P - (A + E))*t[i]) - 1)
            q_t[i] = log(1 - (kf*L*q₀)/(Ω₀*(P - (A + E))) * (exp((P - (A + E))*t[i]) - 1)) / (kf* log(t[i]/Tb))
        end
    end
    return q_t
end

function Ω_t(N₀, P, A, E, kf, Ω₀, dt, t)
    A_t = zeros(length(t))
    for i in eachindex(t)
        A_t[i] = Ω₀ - N_analytic(N₀, P, A, E, t[i]) .* kf .* dt
    end
    return A_t
end

