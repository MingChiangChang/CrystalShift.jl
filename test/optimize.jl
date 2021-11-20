# TODO A lot more test should be written
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
mean_θ = [1., 1., .2]
std_θ = [.025, 10., 1.]
test_path = "data/Ta-Sn-O/sticks.csv"
f = open(test_path, "r")
s = split(read(f, String), "#\n") # Windows: #\r\n ...
#s = Vector{CrystalPhase}()
cs = CrystalPhase.(String.(s[1:end-1]))
# pyro = CrystalPhase(String(s[1]))
# delta = CrystalPhase(String(s[2]))
# push!(cs, pyro)
# push!(cs, delta)
x = collect(8:.1:60)

function test_reconstruct(cp::CrystalPhase, x::AbstractVector, plt=false)
    y = synthesize_data(cp, x)
    c = optimize!(cp, x, y, std_noise, mean_θ, std_θ; maxiter=32, regularization=true)

    if plt
        p = plot(x, c.(x), label="Reconstructed")
        plot!(x, y, label="Original")
        display(p)
    end
    norm(c.(x).-y)
end

function synthesize_data(cp::CrystalPhase, x::AbstractVector)
    params = get_free_params(cp)
    scaling = (0.02.*rand(size(params)[1]).-0.01).+1
    @. params = params*scaling
    params = [params..., 1., 0.2]
    r = reconstruct!(cp, params, x)
    r/max(r...)
end

# test fails when the lattice parameters shift too much
@testset "test" begin
    for (idx, cp) in enumerate(cs)
        println(idx)
        @test test_reconstruct(cp, x, true) < 0.1
    end
end

# test_reconstruct(cs[8], x, true)# 


# TODO use phases with less symmetry
# y = ( reconstruct!(pyro, [10., 0.2, .3], x)
#     + reconstruct!(delta, [5.4, 0.6, .1], x) )
# y = ( reconstruct!(pyro, [17.11, 4.87, 5.5, pi/2, 1., .2], x)
#      + reconstruct!(delta, [4.1, 1., .1], x) )
# y /= max(y...)
# #y = reconstruct!(cs, [10., 0.2, .3, 5.4, 0.6, .1], x)
# # plot(x,cs(x))
# #plot(x, y)

# # Compare old phase mapping and current ones
# #
# res = optimize!(cs, x, y, std_noise, mean_θ, std_θ; maxiter=32, regularization=false)
# println("Res:$(norm(res.(x)-y))")
# plot(x, res.(x), label="reconstruct") # This should give the fitted spectrum
# plot!(x, y, label="real")
