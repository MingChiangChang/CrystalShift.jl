# TODO A lot more test should be written
include("../src/util.jl")
include("../src/peak.jl")
include("../src/Crystal.jl")
include("../src/crystalphase.jl")
include("../src/optimize.jl")

using Plots
using PhaseMapping

test_path = "data/sticks.csv"
f = open(test_path, "r")
s = split(read(f, String), "#\n") # Windows: #\r\n ...
cs = Vector{CrystalPhase}()

pyro = CrystalPhase(String(s[1]))
delta = CrystalPhase(String(s[2]))
push!(cs, pyro)
push!(cs, delta)
x = collect(8:.1:60)

# TODO use phases with less symmetry
y = ( reconstruct!(pyro, [10., 0.2, .3], x)
    + reconstruct!(delta, [5.4, 0.6, .1], x) )
y /= max(y...)
#y = reconstruct!(cs, [10., 0.2, .3, 5.4, 0.6, .1], x)
# plot(x,cs(x))
#plot(x, y)

# Compare old phase mapping and current ones
#
cs = optimize!(cs, x, y; regularization=true)
println("Res:$(norm(cs(x)-y))")
plot(x, cs(x)) # This should give the fitted spectrum
plot!(x, y)
