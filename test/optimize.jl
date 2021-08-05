include("../src/util.jl")
include("../src/peak.jl")
include("../src/crystal.jl")
include("../src/crystalphase.jl")
#include("../src/optimize.jl")

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

# TODO use phases with less symmetry
y = ( reconstruct!(pyro, [10., 0.2, .3], x)
    + reconstruct!(delta, [5.4, 0.6, .1], x) )
y = reconstruct!(cs, [10., 0.2, .3, 5.4, 0.6, .1], x)
#plot(x,cs(x))
plot(x, y)

# Compare old phase mapping and current ones

# optimize!(cs, x, y)
# plot!(x, cs(x)) # This should give the fitted spectrum
