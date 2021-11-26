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

function test_optimize(cp::CrystalPhase, x::AbstractVector, plt=false)
    y = synthesize_data(cp, x)
    c = optimize!(cp, x, y, std_noise, mean_θ, std_θ; maxiter=32, regularization=true)

    if plt
        p = plot(x, c.(x), label="Reconstructed")
        plot!(x, y, label="Original")
        display(p)
    end
    println(c[1].cl)
    println(c[1].origin_cl)
    norm(c.(x).-y)
end

function synthesize_data(cp::CrystalPhase, x::AbstractVector)
    params = get_free_params(cp)
    scaling = (0.05.*rand(size(params, 1)).-0.01).+1
    @. params = params*scaling
    params = [params..., 1., 0.2]
    r = reconstruct!(cp, params, x)
    r/max(r...)
end

function synthesize_multiphase_data(cps::AbstractVector{<:CrystalPhase},
                                   x::AbstractVector)
    r = zero(x)
    full_params = Float64[]
    for cp in cps
        params = get_free_params(cp)
        scaling = (0.005.*rand(size(params, 1)).-0.01).+1
        @. params = params*scaling
        params = vcat(params, 0.5.+3rand(1), 0.1.+0.5(rand(1)))
        full_params = vcat(full_params, params)
        # println(full_params)
    end
    r = reconstruct!(cps, full_params, x)
    r/max(r...)
end

function test_multiphase_optimize(cps::AbstractVector{<:CrystalPhase},
                                   x::AbstractVector, num_phase::Int, plt=false)
    phase = rand(1:size(cps, 1), num_phase)
    y = synthesize_multiphase_data(cps[phase], x)
    c = optimize!(cps[phase], x, y, std_noise, mean_θ, std_θ; maxiter=32, regularization=true)

    if plt
        p = plot(x, c.(x), label="Reconstructed")
        plot!(x, y, label="Original")
        display(p)
    end
    norm(c.(x).-y)
end

# test fails when the lattice parameters shift too much
@testset "Single phase with shift test" begin
    for (idx, cp) in enumerate(cs)
        @test test_optimize(cp, x, true) < 0.1
    end
end

@testset "Multiple phases with shift test" begin
    for _ in 1:10
        @test test_multiphase_optimize(cs, x, 2, true) < 0.1
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
