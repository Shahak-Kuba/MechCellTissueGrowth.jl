using KubaPhD
using Plots
include("../PlottingFncsPDE.jl")
include("FVM_ContinuumSolver.jl")
include("FVM_SolverFncs.jl")

# User input variables
D = 10
kf = 100.0
P = 0.1
A = 0.00
E = 0.05
ρ₀ = 0.05
Tmax = 10
r₀ = 159.15494309189535
btype = "circle"
growth_dir = "inward"

# running simulation
θ,r,ρ,κ,σ,η = FVM_SolveContinuum_Polar(D,kf,P,A,E,ρ₀,Tmax,r₀,btype,growth_dir);

ρ_mean = zeros(size(ρ,1))
for row in eachindex(ρ[:,1])
    ρ_mean[row] = sum(ρ[row,:])/size(ρ,2)
end

# plotting solution
cmap = :jet
p = plotContinuumResults_Polar(θ, r, ρ, cmap, (0,0.2), 180)