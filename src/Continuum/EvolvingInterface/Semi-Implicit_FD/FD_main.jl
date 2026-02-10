using KubaPhD
using LinearAlgebra
include("FD_ContinuumSolvers.jl")
include("../PlottingFncsPDE.jl")
include("FD_SolverFncs.jl")

# simulation parameters
D = 0.15;
kf = 0.01;
A = 0.00;
ρ₀ = 17.18247918322778;
growth_dir = "inward"
Tmax = 30.0
Xmax = 2π

x,h,ρ = FD_SolveContinuumLim_Cartesian(D,kf,A,ρ₀,growth_dir,Tmax,Xmax);

cmap = :jet
f1 = plotContinuumResults_Cartesian(x, h, ρ, D, kf, cmap)


myR = 100ones(120)
T = eltype(myR)
D = 1000;
kf = 20.0;
A = 0.00;
ρ₀ = 0.05;
growth_dir = "inward"
Tmax = 22

θ,R,ρ = FD_SolveContinuumLim_Polar(T,D,kf,A,ρ₀,Tmax,growth_dir,myR);


# simulation parameters
using ForwardDiff
using BackwardDiff

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
out = similar(x)
cfg = ForwardDiff.GradientConfig(r_to_output,x,ForwardDiff.Chunk{20}());

@time ForwardDiff.gradient!(out,r_to_output,x,cfg)

@time ForwardDiff.gradient(r_to_output,x)

# plotContinuumResults_Polar(θ, R, ρ, cmap, D, kf)

