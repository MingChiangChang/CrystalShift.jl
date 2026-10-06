module TestLinearPeakMod
using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate!, evaluate_residual!, get_free_params, get_param_nums
using CrystalShift: PeakModCP, LinearPeakMod, FixedPseudoVoigt, PhaseModel, BackgroundModel, full_optimize!
using CrystalShift: get_PeakModCP, get_newton_objective_func, OptimizationSettings, Simple, EM, DEFAULT_TOL
using CrystalShift: copy_peaks, change_peak_int!
using LinearAlgebra
using CovarianceFunctions: EQ
using Random
using Test

Random.seed!(1)
std_noise = 0.1
mean_θ = [1., .5, .2]
std_θ = [.05, 2., 1.]

test_path = "../data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]))
x = collect(8:.1:60)

heights(IMs) = reduce(vcat, [get_free_params(IM) for IM in IMs])
# Peaks outside of x, or with only a tail inside it, have (nearly) zero basis columns and
# are not determined by the data, only by the prior. Compare heights only for the others.
function identifiable(IMs)
    norms = reduce(vcat, [norm.(eachcol(IM.basis)) for IM in IMs])
    norms .> 1e-2 * maximum(norms)
end
# objective of the analytic route evaluated directly on the pattern
function objective(IMs, y, w, mean, std)
    P = LinearPeakMod(IMs, y, mean, std)
    B = reduce(hcat, [IM.basis for IM in IMs])
    c = sum(IM.const_basis for IM in IMs)
    sum(abs2, y .- c .- B * w) / (2std_noise^2) + sum((log.(w) .- P.mean_log_θ).^2 ./ (2P.std_θ.^2))
end
# synthetic pattern y = Σ (B w_true + c) from PeakModCPs with known height factors
truth(IMs, ws) = PeakModCP.(IMs, ws)
pattern(IMs) = evaluate!(zero(x), PhaseModel(IMs), x)

# objective of the autodiff (BFGS) route, used as the reference
function bfgs_objective(IMs, y, M, mean, std)
    opt = OptimizationSettings{Float64}(std_noise, mean, std, 32, true, bfgs, "LS",
                                        Simple, 1, 1., false, DEFAULT_TOL)
    get_newton_objective_func(PhaseModel(IMs), x, y, opt)(log.(heights(M)))
end

@testset "evaluate_residual! with θ" begin
    IM = PeakModCP(cs[1], x, 10)
    θ = 0.5 .+ rand(10)
    r = ones(length(x))
    evaluate_residual!(IM, θ, x, r)
    @test r ≈ ones(length(x)) .- evaluate!(zero(x), PeakModCP(IM, θ), x)
end

@testset "single phase recovery" begin
    IMs = [PeakModCP(cs[1], x, 10)]
    w_true = 0.5 .+ rand(10)
    y = pattern(truth(IMs, [w_true]))
    P = LinearPeakMod(IMs, y, [1.], [1e3]) # weak prior
    @test get_param_nums(P) == 10
    M = optimize!(P, std_noise; maxiter = 100)
    @test length(M) == 1
    id = identifiable(IMs)
    @test heights(M)[id] ≈ w_true[id] rtol = 1e-3
    @test objective(IMs, y, heights(M), [1.], [1e3]) <= objective(IMs, y, w_true, [1.], [1e3])
    @test norm(pattern(M) - y) < 1e-3
end

@testset "agrees with BFGS route" begin
    for i in (1, 2, 4)
        IMs = [PeakModCP(cs[i], x, 32)]
        n = get_param_nums(IMs[1])
        y = pattern(truth(IMs, [0.5 .+ rand(n)]))
        y .+= 0.01 .* randn(length(x))
        M_new = optimize!(LinearPeakMod(IMs, y, [1.], [.5]), std_noise; maxiter = 32)
        M_ref = optimize!(deepcopy(IMs), x, y, std_noise, [1.], [.5];
                          method = bfgs, objective = "LS", maxiter = 2000, regularization = true)
        F_new = bfgs_objective(IMs, y, M_new, [1.], [.5])
        F_ref = bfgs_objective(IMs, y, M_ref, [1.], [.5])
        @test F_new <= F_ref * (1 + 1e-6)
        @test heights(M_new) ≈ heights(M_ref) rtol = 1e-2
    end
end

@testset "multiple phases" begin
    IMs = [PeakModCP(cs[1], x, 10), PeakModCP(cs[2], x, 32)]
    n2 = get_param_nums(IMs[2])
    ws = [0.5 .+ rand(10), 0.5 .+ rand(n2)]
    y = pattern(truth(IMs, ws))
    for (mean, std) in (([1.], [1e3]), ([1., 1.], [1e3, 1e3])) # shared and per-phase priors
        P = LinearPeakMod(IMs, y, mean, std)
        @test get_param_nums(P) == 10 + n2
        M = optimize!(P, std_noise; maxiter = 100)
        @test length(M) == 2
        @test get_param_nums.(M) == [10, n2]
        id = identifiable(IMs)
        @test heights(M)[id] ≈ vcat(ws...)[id] rtol = 1e-3
        @test objective(IMs, y, heights(M), mean, std) <= objective(IMs, y, vcat(ws...), mean, std)
    end
end

@testset "background in constant term" begin
    bg = BackgroundModel(x, EQ(), 10, rank_tol = 1e-3)
    bg_pattern = evaluate!(zero(x), bg, x)
    IMs = get_PeakModCP(PhaseModel(cs[1:1], nothing, bg), x, 32)
    @test IMs[1].const_basis ≈ PeakModCP(cs[1], x, 32).const_basis .+ bg_pattern
    n = get_param_nums(IMs[1])
    w_true = 0.5 .+ rand(n)
    y = pattern(truth(IMs, [w_true]))
    M = optimize!(LinearPeakMod(IMs, y, [1.], [1e3]), std_noise; maxiter = 100)
    id = identifiable(IMs)
    @test heights(M)[id] ≈ w_true[id] rtol = 1e-3
end

@testset "full_optimize! routing" begin
    pmcp = PeakModCP(CrystalPhase(String(s[1]), 0.1, FixedPseudoVoigt(0.01)), x, 10)
    pmcp = PeakModCP(pmcp, (0.5 .+ 0.5 .* rand(10)) .* get_free_params(pmcp))
    y = evaluate!(zero(x), pmcp, x)
    kw = (method = LM, regularization = true, loop_num = 4, peak_shift_iter = 32,
          mod_peak_num = 32, peak_mod_mean = [1.], peak_mod_std = [.5], peak_mod_iter = 32)
    c_new = full_optimize!(PhaseModel(cs[1:1]), x, y, std_noise, mean_θ, std_θ;
                           objective = "LS", analytic_peak_mod = true, kw...)
    c_old = full_optimize!(PhaseModel(cs[1:1]), x, y, std_noise, mean_θ, std_θ;
                           objective = "LS", analytic_peak_mod = false, kw...)
    r_new = norm(y - evaluate!(zero(x), c_new, x))
    r_old = norm(y - evaluate!(zero(x), c_old, x))
    @test r_new < 0.2
    @test isapprox(r_new, r_old, rtol = 1e-2)

    # modes the analytic route does not cover fall back to the BFGS route
    c_em = full_optimize!(PhaseModel(cs[1:1]), x, y, std_noise, mean_θ, std_θ;
                          objective = "LS", optimize_mode = EM, kw..., loop_num = 1)
    @test c_em isa PhaseModel
    c_kl = full_optimize!(PhaseModel(cs[1:1]), x, y, std_noise, mean_θ, std_θ;
                          objective = "KL", kw..., method = bfgs, loop_num = 1)
    @test c_kl isa PhaseModel
end


@testset "full_optimize! peak-height prior is relative to the reference" begin
    # first 10 peaks twice as high as the reference; a tight prior keeps the fit closer
    cp = copy_peaks(cs[1])
    change_peak_int!(cp, fill(2., 10))
    y = cp.(x)
    ref = [p.I for p in cs[1].peaks]
    # mean log deviation of the fitted from the reference heights of the doubled peaks
    function deviation(loop_num, analytic)
        c = full_optimize!(PhaseModel(cs[1:1]), x, y, std_noise, mean_θ, std_θ;
                           method = LM, loop_num = loop_num, peak_mod_std = [.05],
                           analytic_peak_mod = analytic)
        sum(log.([p.I for p in c.CPs[1].peaks[1:10]] ./ ref[1:10])) / 10
    end
    for analytic in (true, false)
        d1, d8 = deviation(1, analytic), deviation(8, analytic)
        @test 0 < d1 < log(2) / 2 # pulled towards 2, held back by the prior
        # more loops must not loosen the prior (before: per-loop factors compounded, 0.10 -> 0.24)
        @test d8 <= d1 + 0.01
    end
end

end
