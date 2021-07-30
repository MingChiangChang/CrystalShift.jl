include("../src/util.jl")
include("../src/peak.jl")
include("../src/crystal.jl")
include("../src/crystalphase.jl")

using Plots
using PhaseMapping

pyro = CrystalPhase("data/Bi2Ti2O7.csv")

x = collect(8:.1:45)

plot(collect(8:.1:45),pyro(collect(8:.1:45)))
