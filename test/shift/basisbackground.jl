module testbasisbackground
using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, evaluate!, get_free_params
using CrystalShift: get_free_lattice_params, PhaseModel, PseudoVoigt
using CrystalShift: _prior, get_param_nums, lm_prior!, reconstruct_BG!
using CrystalShift: BasisBackground, PolynomialBackground, ChebyshevBackground
using CrystalShift: CosineBackground, InterpolationBackground

using LinearAlgebra
using Random: rand
using Test

x = collect(8:.1:60)

@testset "basis construction" begin
    t = @. 2 * (x - 8) / 52 - 1

    P = PolynomialBackground(x, 3)
    @test get_param_nums(P) == 4
    @test P.basis ≈ hcat(ones(length(x)), t, t.^2, t.^3)
    @test P.λ == [0., 1., 1., 1.] # constant term unregularized by default
    Pinv = PolynomialBackground(x, 2; inverse = true, λ = 2.)
    @test get_param_nums(Pinv) == 4
    @test Pinv.basis[:, end] ≈ 8 ./ x
    @test Pinv.λ == [0., 2., 2., 2.]

    C = ChebyshevBackground(x, 5)
    @test get_param_nums(C) == 6
    for k in 0:5 # Tₖ(cos φ) = cos(kφ)
        @test C.basis[:, k+1] ≈ cos.(k .* acos.(clamp.(t, -1, 1)))
    end

    F = CosineBackground(x, 4)
    @test F.basis[:, 1] ≈ ones(length(x))
    @test F.basis[1, :] ≈ ones(5)
    @test F.basis[end, :] ≈ [(-1)^i for i in 0:4]

    I = InterpolationBackground(x, [10., 30., 50.])
    @test get_param_nums(I) == 3
    @test all(sum(I.basis, dims = 2) .≈ 1) # partition of unity
    I.c .= [1., 3., 2.]
    y = evaluate!(zero(x), I, x)
    @test y[x .<= 10] ≈ fill(1., count(x .<= 10))
    @test y[x .>= 50] ≈ fill(2., count(x .>= 50))
    @test y[findfirst(≈(20.), x)] ≈ 2.
    @test y[findfirst(≈(40.), x)] ≈ 2.5
    @test get_param_nums(InterpolationBackground(x, 7)) == 7

    @test_throws ArgumentError InterpolationBackground(x, [10., 10.])
    @test_throws ArgumentError PolynomialBackground(x .- 20, 2; inverse = true)
    @test_throws DimensionMismatch evaluate!(zeros(3), C, x)
end

@testset "parameters and prior" begin
    B = ChebyshevBackground(x, 3; λ = 2., λ0 = .5)
    θ = [1., 2., 3., 4., 5., 6.]
    rest, new_B = reconstruct_BG!(θ, B)
    @test rest == [5., 6.]
    @test get_free_params(new_B) == [1., 2., 3., 4.]
    @test evaluate!(zero(x), new_B, x) ≈ B.basis * [1., 2., 3., 4.]
    p = zeros(4)
    lm_prior!(p, B, [1., 1., 1., 1.])
    @test p == [.5, 2., 2., 2.]
    @test _prior(B, [1., 1., 1., 1.]) ≈ .25 + 12
end

# joint refinement of a phase together with each background model
verbose = false
maxiter = 2000
std_noise = 1e-3
mean_θ = [1., 1, .2]
std_θ = [.02, 1., 1.]

test_path = "../data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")
if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end
close(f)

cs = @. CrystalPhase(String(s[1:end-1]), (0.1, ), (PseudoVoigt(0.5), ))

function synthesize_data(cp::CrystalPhase, x::AbstractVector)
    params = get_free_lattice_params(cp)
    interval_size = 0.025
    scaling = (interval_size.*rand(size(params, 1),) .- interval_size/2) .+ 1
    @. params = params*scaling
    params = [params..., 1., 0.2, 0.5]
    r = evaluate(cp, params, x)
    normalization = maximum(r)
    params[end-1] /= normalization
    return r/normalization, params
end

data, _ = synthesize_data(cs[1], x)
smooth_bg = @. 0.3 + 0.2 * exp(-(x - 8) / 15) + 0.1 * (x / 60)^2 # air scatter + slow rise
noisy_data = data .+ smooth_bg .+ 0.05 .* rand(length(x))
noisy_data ./= maximum(noisy_data)

@testset "joint refinement with $(name)" for (name, BG) in [
        ("polynomial", PolynomialBackground(x, 4; inverse = true)),
        ("chebyshev", ChebyshevBackground(x, 6)),
        ("cosine", CosineBackground(x, 6)),
        ("interpolation", InterpolationBackground(x, 8)),
    ]
    pm = PhaseModel([cs[1]], nothing, BG)
    c = optimize!(pm, x, noisy_data, std_noise, mean_θ, std_θ,
                  objective = "LS", method = LM, maxiter = maxiter, optimize_mode = Simple,
                  regularization = true, verbose = verbose)
    # noise floor is ~1e-4, a phase-only fit gives ~5e-3
    @test norm(noisy_data .- evaluate!(zero(x), c, x))^2 / length(x) < 1e-3
    @test any(!iszero, get_free_params(c.background))
end

# lm_optimize! fills the Jacobian columns of linear backgrounds from their basis
using CrystalShift: get_lm_objective_func, background_jacobian, SplitJacobianLM, is_linear
using CrystalShift: OptimizationSettings, FixedBackground, BackgroundModel, Wildcard, Lorentz
using OptimizationAlgorithms: update_jacobian!
using CovarianceFunctions: EQ
using ForwardDiff

@testset "constant background Jacobian matches ForwardDiff: $(name)" for (name, BG) in [
        ("basis", ChebyshevBackground(x, 6)),
        ("fixed", FixedBackground(smooth_bg, 1., 3.)),
        ("kernel", BackgroundModel(x, EQ(), 10)),
    ], regularization in (true, false), wildcard in (false, true)
    @test is_linear(BG)
    W = wildcard ? [Wildcard([20.], [.5], [2.], "Amorphous", Lorentz(), [1., 1., .5])] : nothing
    pm = PhaseModel([cs[1]], W, BG)
    opt_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 10, regularization)
    y_uncer = 0.01 .* rand(length(x))
    f = get_lm_objective_func(pm, x, noisy_data, y_uncer, opt_stn)
    nbg = get_param_nums(BG)
    np = get_param_nums(pm) - nbg
    θ = [log.(get_free_params(pm)[1:np]) .+ 0.01 .* randn(np); randn(nbg)]
    r = zeros(regularization ? length(x) + length(θ) : length(x))

    J_ad = ForwardDiff.jacobian(f, copy(r), θ)
    LM = SplitJacobianLM(f, θ, r, np, background_jacobian(BG, x, noisy_data, y_uncer, length(r), opt_stn))
    v, J = update_jacobian!(LM, θ)
    @test J ≈ J_ad rtol = 1e-10
    @test v ≈ f(copy(r), θ)
end

@test !is_linear(nothing)

end
