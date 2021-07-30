include("../src/util.jl")
include("../src/peak.jl")
include("../src/crystal.jl")
include("../src/crystalphase.jl")

using Plots
using PhaseMapping

test_path = "data/sticks.csv"
f = open(test_path, "r")
s = split(read(f, String), '#')
println(s)
cs = Vector{CrystalPhase}()


pyro = CrystalPhase(String(s[1]))

x = collect(8:.1:60)

plot(collect(8:.1:60),pyro(collect(8:.1:60)))
