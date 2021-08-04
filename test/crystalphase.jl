include("../src/util.jl")
include("../src/peak.jl")
include("../src/crystal.jl")
include("../src/crystalphase.jl")

using Plots
using PhaseMapping

test_path = "data/sticks.csv"
f = open(test_path, "r")
s = split(read(f, String), "#\r\n")
println(s[2])
cs = Vector{CrystalPhase}()


pyro = CrystalPhase(String(s[1]))
delta = CrystalPhase(String(s[2]))
push!(cs, pyro)
push!(cs, delta)
x = collect(8:.1:60)

plot(collect(8:.1:60),cs(collect(8:.1:60)))
plot!(collect(8:.1:60), reconstruct!(pyro, [10., 1., .1], collect(8:.1:60)))
