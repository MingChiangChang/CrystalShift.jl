# Robustness + speed of optimization backends on the randomized scenarios from test/.
#
#   julia --project=benchmark benchmark/scenarios.jl [trials-multiplier] [--set packages|new] [--filter a,b]
#
# --set packages (default): the off-the-shelf solvers from solvers.jl
# --set new:                the strategies from methods.jl (multi-start, continuation, VarPro, …)
# --filter a,b:             only solvers whose name contains one of the substrings (the
#                           first solver, the reference for "≤ cur", is always kept)
#
# Each scenario reproduces the data generation, priors, maxiter and pass criterion of one
# test file, but with a fixed seed and many trials. For every solver we report:
#   pass    – fraction of trials meeting that test's own pass criterion
#   lat ok  – fraction of trials whose fitted lattice parameters are all within 0.1% of the
#             true (generating) values; lattice error is max_i |fit_i/true_i - 1| over every
#             free lattice parameter of every phase (synthetic scenarios only)
#   ≤ cur   – fraction of trials whose final cost is no worse than the current LM's (+1%)
#   time    – median / total wall time over trials (second, fully compiled pass)
using Pkg
Pkg.develop(path = joinpath(@__DIR__, ".."); io = devnull)
Pkg.instantiate(; io = devnull)

include(joinpath(@__DIR__, "solvers.jl"))
include(joinpath(@__DIR__, "methods.jl"))
using CrystalShift: BackgroundModel, FixedBackground, Wildcard, FixedPseudoVoigt, PseudoVoigt, Lorentz,
                    evaluate, get_free_lattice_params
using CovarianceFunctions: EQ
using NPZ
using Random
using Statistics
using Printf

const SET = (i = findfirst(==("--set"), ARGS)) === nothing ? "packages" : ARGS[i+1]
const FILTER = (i = findfirst(==("--filter"), ARGS)) === nothing ? nothing : split(ARGS[i+1], ",")
const MULT = something(tryparse.(Int, ARGS)..., 1)
selected(solvers) = FILTER === nothing ? solvers :
    [s for (i, s) in enumerate(solvers) if i == 1 || any(occursin(f, first(s)) for f in FILTER)]
const DATA = joinpath(@__DIR__, "..", "data")
const X = collect(8:0.1:60)

function read_phases(profile)
    s = split(replace(read(joinpath(DATA, "Ta-Sn-O", "sticks.csv"), String), "\r\n" => "\n"), "#\n")
    [CrystalPhase(String(si), 0.1, profile) for si in s if !isempty(strip(si))]
end
const CS_FPV = read_phases(FixedPseudoVoigt(0.5))
const CS_PV = read_phases(PseudoVoigt(0.5))

# --- data generators, copied from the test files ------------------------------------------

# Each returns (y, truth) where truth holds the generating lattice parameters per phase.

# test/shift/optimize.jl, test/shift/background.jl, test/shift/fixedbackground.jl
function synthesize_data(rng, cp, x; profile_param = rand(rng))
    params = get_free_lattice_params(cp)
    interval_size = 0.025
    scaling = (interval_size .* rand(rng, length(params)) .- interval_size / 2) .+ 1
    lattice = params .* scaling
    r = evaluate(cp, [lattice..., 1.0, 0.2, profile_param], x)
    r / maximum(r), [lattice]
end

# test/shift/optimize.jl
function synthesize_multiphase_data(rng, cps, x)
    θ, truth = Float64[], Vector{Float64}[]
    for cp in cps
        params = get_free_lattice_params(cp)
        scaling = (0.01 .* rand(rng, length(params)) .- 0.005) .+ 1
        push!(truth, params .* scaling)
        append!(θ, truth[end], 0.5 + 3rand(rng), 0.1 + 0.1rand(rng))
    end
    r = evaluate(cps, θ, x)
    r / maximum(r), truth
end

# Fitted lattice parameters per phase, from a log-space θ
function fitted_lattices(p::Problem, log_θ)
    out, off = Vector{Float64}[], 0
    for cp in p.pm.CPs
        push!(out, exp.(log_θ[off+1:off+cp.cl.free_param]))
        off += get_param_nums(cp)
    end
    out
end

# max relative lattice error over all phases; phases of the same structure are
# interchangeable, so take the best matching between fitted and true phases
function lattice_error(p::Problem, log_θ, truth)
    fit = fitted_lattices(p, log_θ)
    ids = [cp.id for cp in p.pm.CPs]
    perms = length(fit) == 2 && ids[1] == ids[2] ? ([1, 2], [2, 1]) : (collect(eachindex(fit)),)
    minimum(perms) do perm
        maximum(maximum(abs.(fit[perm[k]] ./ truth[k] .- 1)) for k in eachindex(truth))
    end
end

# --- scenarios: (name, maker(rng, i) -> (Problem, truth), pass(problem, model) -> Bool, …) --

model_resnorm(p, m) = norm(evaluate!(zero(p.x), m, p.x) .- p.y)

struct Scenario
    name::String
    make::Function
    pass::Function
    trials::Int
    solvers::Vector
end

single_phase = Scenario("single phase, ±1.25% strain  [test/shift/optimize.jl]",
    function (rng, i)
        cp = CS_FPV[mod1(i, length(CS_FPV))]
        stn = OptimizationSettings{Float64}(0.1, [1.0, 0.5, 0.2], [0.05, 2.0, 1.0], 512, true, LM)
        y, truth = synthesize_data(rng, cp, X)
        Problem(PhaseModel([cp]), X, y, stn), truth
    end,
    (p, m) -> model_resnorm(p, m) < 0.1, 4 * length(CS_FPV), LS_SOLVERS)

two_phases = Scenario("two phases, ±0.5% strain  [test/shift/optimize.jl]",
    function (rng, i)
        cps = CS_FPV[rand(rng, 1:length(CS_FPV), 2)]
        stn = OptimizationSettings{Float64}(0.1, [1.0, 0.5, 0.2], [0.05, 2.0, 1.0], 512, true, LM, "LS",
                                            CrystalShift.Simple, 1, 0.1)
        y, truth = synthesize_multiphase_data(rng, cps, X)
        Problem(PhaseModel(cps), X, y, stn), truth
    end,
    (p, m) -> model_resnorm(p, m) < 0.1, 40, LS_SOLVERS)

with_background = Scenario("single phase + smooth background  [test/shift/background.jl]",
    function (rng, i)
        cp = CS_FPV[mod1(i, length(CS_FPV))]
        y, truth = synthesize_data(rng, cp, X; profile_param = 0.5)
        y = y .+ 0.1 .* (1 .+ sin.(0.2 .* X))
        bg = BackgroundModel(X, EQ(), 10, rank_tol = 1e-4)
        stn = OptimizationSettings{Float64}(1e-3, [1.0, 1, 0.2], [0.02, 1.0, 1.0], 2000, true, LM)
        Problem(PhaseModel([cp], nothing, bg), X, y, stn), truth
    end,
    (p, m) -> model_resnorm(p, m) < 0.1, length(CS_FPV), LS_SOLVERS)

const POLY = let b = @. (X - 2) * (X - 20) * (X - 30) * (X - 50) * (X - 70)
    b .-= minimum(b); b ./= maximum(b) * 2
end
fixed_background = Scenario("single phase + polynomial bg + noise  [test/shift/fixedbackground.jl]",
    function (rng, i)
        cp = CS_PV[mod1(i, length(CS_PV))]
        y, truth = synthesize_data(rng, cp, X; profile_param = 0.5)
        y = y .+ POLY .+ 0.05 .* rand(rng, length(X))
        y ./= maximum(y)
        stn = OptimizationSettings{Float64}(1e-3, [1.0, 1, 0.2], [0.02, 1.0, 1.0], 2000, true, LM)
        Problem(PhaseModel([cp], nothing, FixedBackground(POLY, 1.0, 10.0)), X, y, stn), truth
    end,
    (p, m) -> model_resnorm(p, m)^2 / length(p.x) < 0.01, length(CS_PV), LS_SOLVERS)

# Measured data, deterministic: one trial. The test itself uses bfgs, so include QN solvers.
wildcard = Scenario("wildcard + background, measured data  [test/shift/wildcard.jl]",
    function (rng, i)
        q = npzread(joinpath(DATA, "test_q.npy"))
        y = npzread(joinpath(DATA, "test_int.npy"))
        y ./= maximum(y) * 2
        y ./= maximum(y) * 2  # fit_amorphous normalizes a second time
        w = Wildcard([20.0, 35.0], [1.0, 0.2], [2.0, 3.0], "Amorphous", Lorentz(), [2.0, 2.0, 1.0, 1.0, 0.2, 0.5])
        bg = BackgroundModel(q, EQ(), 20, rank_tol = 1e-3)
        stn = OptimizationSettings{Float64}(1e-2, [1.0, 1.0, 1.0], [1.0, 1.0, 1.0], 512, true, LM)
        Problem(PhaseModel(w, bg), q, y, stn), nothing
    end,
    (p, m) -> model_resnorm(p, m) < 0.3, 1, [LS_SOLVERS; QN_SOLVERS])

const SCENARIOS = [single_phase, two_phases, with_background, fixed_background, wildcard]

# --- runner -----------------------------------------------------------------------------

function run_scenario(sc::Scenario)
    ntrials = sc.trials == 1 ? 1 : sc.trials * MULT
    rng = Xoshiro(1234)
    made = [sc.make(rng, i) for i in 1:ntrials]
    problems, truths = first.(made), last.(made)
    θ0s = initial_log_θ.(problems)
    @printf("\n%s — %d trials, %d params\n", sc.name, ntrials, length(θ0s[1]))
    @printf("  %-40s %6s %6s %6s %10s %10s %6s %12s %16s\n", "solver", "pass", "lat ok", "≤ cur",
            "median", "total", "iters", "median cost", "lat err med/p90")

    ref_cost = nothing
    for (name, solve) in selected(SET == "new" ? NEW_SOLVERS : sc.solvers)
        # Untimed pass first: phases of different crystal systems are different types, and
        # mixed-system models hit new method combinations, so compile on every trial.
        for (p, θ0) in zip(problems, θ0s)
            try solve(p, θ0) catch end
        end
        times, costs, iters, laterr, passed, errors = Float64[], Float64[], Int[], Float64[], 0, 0
        for (p, θ0, truth) in zip(problems, θ0s, truths)
            local θ, it
            t = @elapsed try
                θ, it = solve(p, θ0)
            catch
                θ, it = nothing, missing
            end
            push!(times, t)
            if θ === nothing || any(!isfinite, θ)
                errors += 1; push!(costs, Inf); push!(laterr, Inf); continue
            end
            truth === nothing || push!(laterr, lattice_error(p, θ, truth))
            it === missing || push!(iters, it)
            push!(costs, ls_cost(p, θ))
            passed += sc.pass(p, to_model(p, θ))
        end
        ref_cost === nothing && (ref_cost = costs)
        no_worse = count(costs .<= ref_cost .* 1.01 .+ 1e-12) / ntrials
        latok = isempty(laterr) ? "-" : @sprintf("%5.0f%%", 100count(<(1e-3), laterr) / ntrials)
        latstat = isempty(laterr) ? "-" : @sprintf("%.3f%% / %.3f%%", 100median(laterr), 100quantile(laterr, 0.9))
        @printf("  %-40s %5.0f%% %6s %5.0f%% %10s %10s %6s %12.5g %16s%s\n", name, 100passed / ntrials, latok,
                100no_worse, BenchmarkTools.prettytime(1e9median(times)), BenchmarkTools.prettytime(1e9sum(times)),
                isempty(iters) ? "-" : string(round(Int, median(iters))), median(costs), latstat,
                errors > 0 ? "  ($errors errored)" : "")
    end
end

using BenchmarkTools: BenchmarkTools
println("Julia ", VERSION, ", threads ", Threads.nthreads())
println("Costs are ‖get_lm_objective_func residual‖² for every solver (including quasi-Newton ones).")
foreach(run_scenario, SCENARIOS)
