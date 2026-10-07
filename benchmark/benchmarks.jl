# Benchmark suite for CrystalShift. Defines `SUITE` (BenchmarkTools.BenchmarkGroup),
# so it can be run by `benchmark/run.jl` or by PkgBenchmark.
#
# All inputs are deterministic (fixed seed, fixed lattice perturbation) so that the
# optimizer does the same amount of work from run to run. If a code change alters
# numerics, the iteration count can change too, so check the `checks` printed by run.jl.
using BenchmarkTools
using CrystalShift
using CrystalShift: CrystalPhase, PhaseModel, BackgroundModel, OptimizationSettings,
                    FixedPseudoVoigt, Lorentz, LM, Simple, EM,
                    evaluate, evaluate!, evaluate_residual!, get_free_params,
                    get_free_lattice_params, get_lm_objective_func, extend_priors
using CovarianceFunctions: EQ
using ForwardDiff
using LinearAlgebra
using Random

const DATA = joinpath(@__DIR__, "..", "data")

function read_phases(path, profile)
    s = split(replace(read(path, String), "\r\n" => "\n"), "#\n")
    [CrystalPhase(String(si), 0.1, profile) for si in s if !isempty(strip(si))]
end

# Synthetic pattern from `cps` with every free lattice parameter scaled by `strain`.
function synthesize(cps, x; strain = 1.005, act = 1.0, σ = 0.2)
    θ = Float64[]
    for cp in cps
        append!(θ, get_free_lattice_params(cp) .* strain, act, σ)
    end
    y = evaluate(cps, θ, x)
    y ./ maximum(y)
end

Random.seed!(0)

# Ta-Sn-O phase library, same as test/optimize.jl
const TASNO = read_phases(joinpath(DATA, "Ta-Sn-O", "sticks.csv"), FixedPseudoVoigt(0.5))
const X = collect(8:0.1:60)
const ONE = [TASNO[10]]
const THREE = TASNO[[2, 5, 10]]
const Y1 = synthesize(ONE, X)
const Y3 = synthesize(THREE, X)

const STD_NOISE = 0.1
const MEAN_θ = [1.0, 0.5, 0.2]
const STD_θ = [0.05, 2.0, 1.0]
const MAXITER = 128

opt_kw(; kw...) = (; method = LM, maxiter = MAXITER, regularization = true, verbose = false, kw...)

# Run the LM objective as the optimizer sees it (log-space θ, regularized residual)
function lm_objective(pm, x, y)
    stn = OptimizationSettings{Float64}(STD_NOISE, MEAN_θ, STD_θ, MAXITER)
    f = get_lm_objective_func(pm, x, y, zero(y), stn)
    log_θ = log.(get_free_params(pm))
    r = zeros(length(y) + length(log_θ))
    f, r, log_θ
end

const SUITE = BenchmarkGroup()

# --- forward model -------------------------------------------------------------
fw = SUITE["forward"] = BenchmarkGroup()
fw["evaluate 1 phase"] = @benchmarkable evaluate!(y, $(ONE[1]), $X) setup = (y = zero($X))
fw["evaluate 3 phases"] = @benchmarkable evaluate!(y, $THREE, $X) setup = (y = zero($X))
let θ = get_free_params(THREE)
    fw["evaluate_residual 3 phases"] = @benchmarkable evaluate_residual!($THREE, $θ, $X, r) setup = (r = copy($Y3))
end

# --- objective + Jacobian (inner loop of LM) ----------------------------------
obj = SUITE["objective"] = BenchmarkGroup()
for (name, cps, y) in (("1 phase", ONE, Y1), ("3 phases", THREE, Y3))
    f, r, log_θ = lm_objective(PhaseModel(cps), X, y)
    obj["residual "*name] = @benchmarkable $f($r, $log_θ)
    obj["jacobian "*name] = @benchmarkable ForwardDiff.jacobian($f, $r, $log_θ)
end

# --- full optimization ---------------------------------------------------------
# Phases are immutable, so optimize! never mutates its inputs; no setup needed.
op = SUITE["optimize"] = BenchmarkGroup()
op["LM 1 phase"] = @benchmarkable optimize!($ONE, $X, $Y1, STD_NOISE, MEAN_θ, STD_θ; opt_kw()...)
op["LM 3 phases"] = @benchmarkable optimize!($THREE, $X, $Y3, STD_NOISE, MEAN_θ, STD_θ; opt_kw()...)
op["LM 3 phases EM"] = @benchmarkable optimize!($THREE, $X, $Y3, STD_NOISE, MEAN_θ, STD_θ;
                                                opt_kw(optimize_mode = EM, em_loop_num = 4)...)
let bg = BackgroundModel(X, EQ(), 5.0), y = Y3 .+ 0.1 .* (1 .+ sin.(0.2 .* X))
    op["LM 3 phases + background"] = @benchmarkable optimize!($(PhaseModel(THREE, nothing, bg)), $X, $y,
                                                              STD_NOISE, MEAN_θ, STD_θ; opt_kw()...)
end

# Cheap correctness fingerprints of the optimize outputs. A perf change that
# moves these is a numerical change, and timings are no longer apples to apples.
function checks()
    resnorm(cps, y) = norm(evaluate!(zero(X), cps, X) .- y)
    pm_bg = PhaseModel(THREE, nothing, BackgroundModel(X, EQ(), 5.0))
    y_bg = Y3 .+ 0.1 .* (1 .+ sin.(0.2 .* X))
    Dict(
        "LM 1 phase resnorm" => resnorm(optimize!(ONE, X, Y1, STD_NOISE, MEAN_θ, STD_θ; opt_kw()...), Y1),
        "LM 3 phases resnorm" => resnorm(optimize!(THREE, X, Y3, STD_NOISE, MEAN_θ, STD_θ; opt_kw()...), Y3),
        "LM 3 phases + background resnorm" => norm(evaluate!(zero(X), optimize!(pm_bg, X, y_bg, STD_NOISE, MEAN_θ, STD_θ; opt_kw()...), X) .- y_bg),
    )
end
