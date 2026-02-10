import MechCellTissueGrowth as MCTG

# Simulation parameters
Seed = 123 # set random seed number for reproducability 

# Domain Properties
Domain = MCTG.DomainProperties_t(R₀=50.0, btype="circle", N = 30, m=4)

# Cell mechanics
CellType1 = MCTG.CellMechProperties_t(kₛ=25, growth_dir="inward", a=20.0, kf=50);
HomCellMech = MCTG.generate_homogeneous_population(CellType1, Domain.N, Domain.m);

# Cell Behaviours
Prolif = MCTG.CellEvent_t(true, 0.1, "Constant");
Death = MCTG.CellEvent_t(true, 0.01, "Constant");
Embed = MCTG.CellEvent_t(true, 0.1, "Constant");
ProlifEmbed = MCTG.CellEvent_t(); # simultaneously proliferate and embed

# Time parameters
SimTime = MCTG.SimTime_t(Tmax=5, δt = 0.001, event_δt = 0.001)

# Run simulation
sol, embed_pos, embed_count, embed_rates = MCTG.GrowthSimulation(Domain, HomCellMech, SimTime, Prolif, Death, Embed, ProlifEmbed, Seed, 100);