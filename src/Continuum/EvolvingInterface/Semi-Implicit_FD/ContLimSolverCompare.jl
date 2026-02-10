using KubaPhD
using LinearAlgebra
using ForwardDiff
include("FD_ContinuumSolvers.jl")
include("../PlottingFncsPDE.jl")
include("FD_SolverFncs.jl")
include("FD_SolverOptimized.jl")


myR = 100ones(100)
T = eltype(myR)
D = 1000;
kf = 20.0;
A = 0.00;
ρ₀ = 0.05;
growth_dir = "inward"
Tmax = 22

@time θ_1,R_1,ρ_1 = FD_SolveContinuumLim_Polar(T,D,kf,A,ρ₀,Tmax,growth_dir,myR);
@time θ_2,R_2,ρ_2 = FD_SolveContinuumLim_Polar_Optimized(T,D,kf,A,ρ₀,Tmax,growth_dir,myR);

# check difference
θ_1 - θ_2
R_1 - R_2
ρ_1 - ρ_2


# Check gradient form dual numbers

function r_to_output(myR)
    T = eltype(myR)
    D = 1000;
    kf = 20.0;
    A = 0.00;
    ρ₀ = 0.05;
    growth_dir = "inward"
    Tmax = 22

    θ,R,ρ = FD_SolveContinuumLim_Polar(T,D,kf,A,ρ₀,Tmax,growth_dir,myR);
    ρ[end]
end

x = 100ones(100)
∇_1 = similar(x)
cfg = ForwardDiff.GradientConfig(r_to_output,x,ForwardDiff.Chunk{20}());
@time ForwardDiff.gradient!(∇_1,r_to_output,x,cfg)

function r_to_output_Opt(myR)
    T = eltype(myR)
    D = 1000;
    kf = 20.0;
    A = 0.00;
    ρ₀ = 0.05;
    growth_dir = "inward"
    Tmax = 22

    θ,R,ρ = FD_SolveContinuumLim_Polar_Optimized(T,D,kf,A,ρ₀,Tmax,growth_dir,myR);
    ρ[end]
end

@time ∇_2 = ForwardDiff.gradient(r_to_output_Opt, x)

∇_1 - ∇_2