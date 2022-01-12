# include("../src/CrystalShift.jl")
# TODO: if this should be in the testing suite, add tests and remove plots and @time macros
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
y = zero(x)

# plot(collect(8:.1:60),cs(collect(8:.1:60)))
# plot!(collect(8:.1:60), reconstruct!(pyro, [10., 1., .1], collect(8:.1:60)))
@time a = reconstruct!(pyro, [10., 1., .1], x)
@time a = reconstruct!(delta, [5.3, 1., .1], x)
@time reconstruct!(pyro, [10., 1., .1], x, y)
@time reconstruct!(delta, [5.3, 1., .1], x, y)
plot(x, y)
