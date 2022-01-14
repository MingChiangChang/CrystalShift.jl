module Testoptimize
using CrystalShift
using CrystalShift: CrystalPhase, optimize!, reconstruct!, get_free_params
using CrystalShift: newton!

using LinearAlgebra
using Random: rand
using Test

verbose = false

# Global
std_noise = .01
mean_θ = [1., 1e-4, .2]
std_θ = [.05, 100, 1.]

test_path = "../data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]))
x = collect(8:.1:60)

function test_optimize(cp::CrystalPhase, x::AbstractVector,
                       method::OptimizationMethods, verbose = verbose)
    y, sol = synthesize_data(cp, x)
    c = optimize!(cp, x, y, std_noise, mean_θ, std_θ;
                  method = method,
                  maxiter = 256, regularization = true)

    if verbose
        println(c[1].cl)
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
    interval_size = 0.01

    for cp in cps
        params = get_free_params(cp)

        scaling = (interval_size.*rand(size(params, 1)).-interval_size/2).+1
        @. params = params*scaling
        params = vcat(params, 0.5.+3rand(1), 0.1.+0.1(rand(1)))
        full_params = vcat(full_params, params)
    end
    r = reconstruct!(cps, full_params, x)

    r/maximum(r)
end

function test_multiphase_optimize(cps::AbstractVector{<:CrystalPhase},
                                   x::AbstractVector, num_phase::Int,
                                   method::OptimizationMethods, verbose = false)
    phase = rand(1:size(cps, 1), num_phase)
    y = synthesize_multiphase_data(cps[phase], x)
    c = optimize!(cps[phase], x, y, std_noise, mean_θ, std_θ;
                  method = method, maxiter = 256,
                  regularization = true, verbose = verbose)

    if verbose
        println(c)
    end
    
    norm(c.(x).-y)
end

# test fails when the lattice parameters shift too much
@testset "Single phase with shift LM test" begin
    for (idx, cp) in enumerate(cs)
        verbose && println(idx)
        @test test_optimize(cp, x, LM, false) < 0.1
    end
end

@testset "Single phase with shift newton test" begin
    for (idx, cp) in enumerate(cs)
        verbose && println(idx)
        c = @test test_optimize(cp, x, Newton, false) < 0.1
    end
end

@testset "Multiple phases with shift test" begin
    for _ in 1:5
        @test test_multiphase_optimize(cs, x, 2, LM, false) < 0.1
    end
end

@testset "Multiple phases with shift newton test" begin
    for _ in 1:5
        @test test_multiphase_optimize(cs, x, 2, Newton, false) < 0.1
    end
end

end # module