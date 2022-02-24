module Testoptimize
using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params
using CrystalShift: newton!, get_free_lattice_params, get_fraction

using LinearAlgebra
using Random: rand
using Test
# using BenchmarkTools

verbose = false
residual_tol = 0.1 # tolerance for residual norm after optimization
maxiter = 512 # appears to be required for phase combinations in particular

# Global
std_noise = 0.1
mean_θ = [1., 1, .2]
std_θ = [.02, 1., 1.]
# newton_lambda = 1e-2 TODO: make this passable to the newton optimization

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
                  maxiter = maxiter, regularization = true)

    if verbose
        println(c[1].cl)
        println("sol: $(sol)")
        println(c[1].origin_cl)
    end

    norm(c.(x).-y)
end

function synthesize_data(cp::CrystalPhase, x::AbstractVector)
    params = get_free_lattice_params(cp)
    interval_size = 0.025
    scaling = (interval_size.*rand(size(params, 1),) .- interval_size/2) .+ 1
    @. params = params*scaling
    params = [params..., 1., 0.2, 0.5]
    r = evaluate(cp, params, x)
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
        params = get_free_lattice_params(cp)

        scaling = (interval_size.*rand(size(params, 1)).-interval_size/2).+1
        @. params = params*scaling
        params = vcat(params, 0.5.+3rand(1), 0.1.+0.1(rand(1)), (rand(1)))
        full_params = vcat(full_params, params)
    end
    r = evaluate(cps, full_params, x)

    r/maximum(r)
end

function test_multiphase_optimize(cps::AbstractVector{<:CrystalPhase},
                                   x::AbstractVector, num_phase::Int,
                                   method::OptimizationMethods, objective::String = "LS", 
                                   verbose = false)
    phase = rand(1:size(cps, 1), num_phase)
    y = synthesize_multiphase_data(cps[phase], x)
    c = optimize!(cps[phase], x, y, std_noise, mean_θ, std_θ;
                  objective = objective, method = method, maxiter = maxiter,
                  regularization = true, verbose = verbose)

    if verbose
        println(c)
    end

    norm(c.(x).-y)
end

@time test_optimize(cs[10], x, LM, verbose)
@time test_optimize(cs[10], x, Newton, false)
@time test_multiphase_optimize(cs, x, 2, LM, "LS", verbose)


# test fails when the lattice parameters shift too much
@testset "Single phase with shift LM test" begin
    correct_counts = 0
    for (idx, cp) in enumerate(cs)
        verbose && println(idx)
        if test_optimize(cp, x, LM, verbose) < 0.1
            correct_counts += 1
        end
    end
    @test correct_counts >= size(cs, 1) - 1
end

@testset "Single phase with shift newton test" begin
    correct_counts = 0
    for (idx, cp) in enumerate(cs)
        verbose && println(idx)
        if  test_optimize(cp, x, Newton, false) < 0.1
            correct_counts += 1
        end
    end
    @test correct_counts >= size(cs, 1) - 1
end

# @testset "Single phase with shift bfgs test" begin
#     correct_counts = 0
#     for (idx, cp) in enumerate(cs)
#         verbose && println(idx)
#         if  test_optimize(cp, x, bfgs, true) < 0.1 # FIXME: Sometimes sin(x) goes to infinite (This only happens when doing test)
#             correct_counts += 1                    # FIXME: dx goes to 0.0 when std_noise goes to 1e-2
#         end
#     end
#     @test correct_counts >= size(cs, 1) - 1
# end

# @testset "Single phase with shift lbfgs test" begin
#     correct_counts = 0
#     for (idx, cp) in enumerate(cs)
#         verbose && println(idx)
#         if  test_optimize(cp, x, l_bfgs, true) < 0.1
#             correct_counts += 1
#         end
#     end
#     @test correct_counts >= size(cs, 1) - 1
# end

@testset "Multiple phases with shift test" begin
    correct_counts = 0
    for _ in 1:5
        if test_multiphase_optimize(cs, x, 2, LM, "LS", verbose) < 0.1
            correct_counts += 1
        end
    end
    @test correct_counts >= 3
end

# @testset "Multiple phases with shift newton KL test" begin # call it a success if 4/5 cases passed
#     correct_counts = 0
#     for _ in 1:5
#         if test_multiphase_optimize(cs, x, 2, Newton, "KL", verbose) < 0.1
#             correct_counts += 1
#         end
#     end
#     @test correct_counts >= 3
# end

@testset "Multiple phases with shift newton least sqaure test" begin # call it a success if 4/5 cases passed
    correct_counts = 0
    for _ in 1:5
        if test_multiphase_optimize(cs, x, 2, Newton, "LS", verbose) < 0.1
            correct_counts += 1
        end
    end
    @test correct_counts >= 3
end

end # module
