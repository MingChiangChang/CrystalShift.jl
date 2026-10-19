module TestDogleg
using CrystalShift
using CrystalShift: CrystalPhase, PhaseModel, OptimizationSettings, dogleg, optimize!, get_free_params,
                    get_free_lattice_params, Wildcard, BackgroundModel, EQ, USE_BLOCK_DOGLEG,
                    DOGLEG_MAX_STRAIN
using CrystalShift: Lorentz, Gauss, PseudoVoigt, FixedPseudoVoigt
using LinearAlgebra
using Test

# dogleg_loop! (block Jacobians) against LeastSquaresOptim's Dogleg on the same problem. With
# pred_tol = 0 it is the same algorithm (skipped retries are exact), so the fits agree to
# round-off; with the default pred_tol it stops a little earlier, at the same fit quality.

f = open("../data/tree_test_sticks.csv", "r")
s = Sys.iswindows() ? split(read(f, String), "#\r\n") : split(read(f, String), "#\n")
close(f)
s = filter(!isempty ∘ strip, s)

x = collect(LinRange(8, 45, 800))
pseudo_noise(n, σ) = σ * sqrt(3) .* (2 .* mod.(sin.(12.9898 .* (1:n)) .* 43758.5453, 1) .- 1)
strained(cp, f) = (θ = copy(get_free_params(cp)); θ[1:cp.cl.free_param] .*= f; CrystalPhase(cp, θ))
opt = OptimizationSettings{Float64}(.01, [1., 1., .2], [.5, .5, 1.], 128, true, dogleg)

# fitted model of a fresh model, with the block loop on or off
function fit(make_model, y; block::Bool)
    try
        USE_BLOCK_DOGLEG[] = block
        optimize!(make_model(), x, y, opt)
    finally
        USE_BLOCK_DOGLEG[] = true
    end
end
resnorm(pm, y) = norm(y .- CrystalShift.evaluate!(zero(x), pm, x))

# log-space parameters after dogleg_optimize!, with the block loop on or off
function dogleg_log_θ(pm, y; block::Bool, pred_tol)
    log_θ = log.(get_free_params(pm))
    try
        USE_BLOCK_DOGLEG[] = block
        CrystalShift.dogleg_optimize!(log_θ, pm, x, y, zero(y), opt; pred_tol)
    finally
        USE_BLOCK_DOGLEG[] = true
    end
end

@testset "Block dogleg matches LeastSquaresOptim dogleg" begin
    for profile in (Lorentz(), Gauss(), PseudoVoigt(0.5), FixedPseudoVoigt(0.5)), ids in ([1], [1, 6])
        cs = CrystalPhase.(String.(s), (0.1,), (profile,))
        y = sum(strained(cs[id], k == 1 ? 1.005 : 0.996).(x) for (k, id) in enumerate(ids))
        y = y ./ maximum(y) .+ pseudo_noise(length(x), 0.01)
        model() = PhaseModel(cs[ids])
        # same algorithm without the predicted-decrease stop: same result to round-off
        lso = dogleg_log_θ(model(), y; block = false, pred_tol = 0)
        @test dogleg_log_θ(model(), y; block = true, pred_tol = 0) ≈ lso rtol = 1e-8
        # default: stops a little earlier at the same fit quality
        a, b = fit(model, y; block = true), fit(model, y; block = false)
        @test resnorm(a, y) ≈ resnorm(b, y) rtol = 1e-6
        @test get_free_params(a) ≈ get_free_params(b) rtol = 1e-4
    end
end

@testset "Block dogleg with a linear background and a wildcard" begin
    cs = CrystalPhase.(String.(s), (0.1,), (PseudoVoigt(0.5),))
    y = strained(cs[1], 1.004).(x) .+ 0.2 .* exp.(-((x .- 25) ./ 15) .^ 2)
    y = y ./ maximum(y) .+ pseudo_noise(length(x), 0.01)
    bg() = PhaseModel(cs[[1]], nothing, BackgroundModel(x, EQ(), 8., 10., rank_tol = 1e-3))
    @test resnorm(fit(bg, y; block = true), y) ≈ resnorm(fit(bg, y; block = false), y) rtol = 1e-6
    w() = PhaseModel(cs[[1]], [Wildcard([25.], [.2], [5.], "Amorphous", Lorentz(), [.5, .5, .5])], nothing)
    @test resnorm(fit(w, y; block = true), y) ≈ resnorm(fit(w, y; block = false), y) rtol = 1e-6
end

@testset "Block dogleg lattice bound" begin
    # data 10% strained: each lattice parameter may move at most DOGLEG_MAX_STRAIN
    cs = CrystalPhase.(String.(s), (0.1,), (Lorentz(),))
    y = strained(cs[6], 1.1).(x); y ./= maximum(y)
    c = optimize!(PhaseModel(cs[[6]]), x, y, opt).CPs[1]
    strain = get_free_lattice_params(c) ./ get_free_lattice_params(cs[6]) .- 1
    @test all(abs.(strain) .<= DOGLEG_MAX_STRAIN + 1e-9)
end

@testset "Background-only model uses LeastSquaresOptim" begin
    # no phase parameters: the block path does not apply, the fallback must still fit
    y = 0.5 .+ 0.2 .* sin.(x ./ 7)
    model = PhaseModel(BackgroundModel(x, EQ(), 8., 10., rank_tol = 1e-3))
    r(θ) = norm(y .- CrystalShift.evaluate!(zero(x), PhaseModel(model, θ), x))
    θ = get_free_params(optimize!(model, x, y, opt))
    @test all(isfinite, θ)
    @test r(θ) < r(get_free_params(model))
end

println("End of dogleg.jl test")
end # module
