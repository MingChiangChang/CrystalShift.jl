# TODO A lot more test should be written
using CrystalShift
using CrystalShift: CrystalPhase, optimize!, reconstruct!, get_free_params
using CrystalShift: newton!

using LinearAlgebra
# using Plots # let's only load this when debugging locally
using Random: rand
using Test

verbose = false

# Global
std_noise = .01
mean_θ = [1., 1., .2]
std_θ = [.05, Inf, .1]

test_path = "../data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
# test_path = "data/Ta-Sn-O/sticks.csv"
# test_path = "CrystalShift.jl/data/Ta-Sn-O/sticks.csv"
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]))
x = collect(8:.1:60)

function test_optimize(cp::CrystalPhase, x::AbstractVector, plt=false, verbose = verbose)
    y, sol = synthesize_data(cp, x)
    c = optimize!(cp, x, y, std_noise, mean_θ, std_θ;
                  method = LM,
                  maxiter = 64, regularization = true)

    if plt
        p = plot(x, c.(x), label="Reconstructed")
        plot!(x, y, label="Original")
        display(p)
    end
    if verbose
        println(c[1].cl)
        println("sol: $(sol)")
        println(c[1].origin_cl)
    end
    norm(c.(x).-y)
end

function test_newton(cp::CrystalPhase, x::AbstractVector, plt=false, verbose = verbose)
    y, sol = synthesize_data(cp, x)
    c = optimize!(cp, x, y, std_noise, mean_θ, std_θ;
                  method = Newton,
                  maxiter = 64, regularization = true, verbose = verbose)

    if plt
        p = plot(x, c.(x), label="Reconstructed")
        plot!(x, y, label="Original")
        display(p)
    end

    if verbose
        println(c)
        println("sol: $(sol)")
        println(c[1].origin_cl)
    end
    norm(c.(x).-y)
end

function synthesize_data(cp::CrystalPhase, x::AbstractVector)
    params = get_free_params(cp)
    interval_size = 0.025
    scaling = (interval_size.*rand(size(params, 1),) .- interval_size/2) .+ 1
    @. params = params*scaling
    params = [params..., 1., 0.2]
    r = reconstruct!(cp, params, x)
    normalization = maximum(r)
    params[end-1] /= normalization
    if verbose
        println("synthesize_data: ", params)
    end
    return r/normalization, params
end

function synthesize_multiphase_data(cps::AbstractVector{<:CrystalPhase},
                                   x::AbstractVector)
    r = zero(x)
    full_params = Float64[]
    interval_size = 0.025
    for cp in cps
        params = get_free_params(cp)

        scaling = (interval_size.*rand(size(params, 1)).-interval_size/2).+1
        @. params = params*scaling
        params = vcat(params, 0.5.+3rand(1), 0.1.+0.5(rand(1)))
        full_params = vcat(full_params, params)
        # println(full_params)
    end
    r = reconstruct!(cps, full_params, x)
    r/maximum(r)
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
@testset "Single phase with shift LM test" begin
    for (idx, cp) in enumerate(cs)
        verbose && println(idx)
        @test test_optimize(cp, x, false) < 0.1
    end
end

@testset "Single phase with shift newton test" begin
    for (idx, cp) in enumerate(cs)
        verbose && println(idx)
        c = @test test_newton(cp, x, false) < 0.1
    end
end

# @testset "Multiple phases with shift test" begin
#     for _ in 1:10
#         @test test_multiphase_optimize(cs, x, 2, true) < 0.1
#     end
# end

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
