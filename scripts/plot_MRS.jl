include("../src/util.jl")
include("../src/peak.jl")
include("../src/crystal.jl")
include("../src/crystalphase.jl")
include("../src/optimize.jl")

using Plots
using PhaseMapping
using Random: rand
using Test

# Global
std_noise = .01
mean_θ = [1., .2]
std_θ = [.025, 1.]

test_path = "data/Ta-Sn-O/sticks.csv"
f = open(test_path, "r")
s = split(read(f, String), "#\r\n") # Windows: #\r\n ...
#s = Vector{CrystalPhase}()
cs = CrystalPhase.(String.(s[1:end-1]))
# pyro = CrystalPhase(String(s[1]))
# delta = CrystalPhase(String(s[2]))
# push!(cs, pyro)
# push!(cs, delta)
x = collect(8:.1:60)

i=3
p_pos = Float64[]
p_int = Float64[]

for p in cs[i].peaks
    push!(p_pos, cs[i].cl(p)*10)
    push!(p_int, p.I)
end

plot(x, CrystalPhase(cs[i], [17., 4.8, 5.5, pi/2, 1.0, 0.3]).(x),
        lw=3, legend=false, xtickfontsize=10,ytickfontsize=10)
plot!(p_pos, p_int, st=:sticks, lw=3)
xlabel!("Q (1/nm)")
ylabel!("a.u.")
