# Head-to-head timing of optimization backends on three fixed, easy problems.
#
#   julia --project=benchmark benchmark/optimizers.jl
#
# Stopping rules differ between packages, so compare time together with the final cost:
# a solver that is faster but stops at a higher cost isn't actually better.
# For success rates on the randomized scenarios from test/, see scenarios.jl.
using Pkg
Pkg.develop(path = joinpath(@__DIR__, ".."); io = devnull)
Pkg.instantiate(; io = devnull)

include(joinpath(@__DIR__, "benchmarks.jl"))  # problem data (ONE, THREE, X, Y1, Y3, priors)
include(joinpath(@__DIR__, "solvers.jl"))
using Printf

stn = OptimizationSettings{Float64}(STD_NOISE, MEAN_θ, STD_θ, MAXITER, true, LM)
const PROBLEMS = [
    "1 phase" => Problem(PhaseModel(ONE), X, Y1, stn),
    "3 phases" => Problem(PhaseModel(THREE), X, Y3, stn),
    "3 phases + background" => Problem(PhaseModel(THREE, nothing, BackgroundModel(X, EQ(), 5.0)), X,
                                       Y3 .+ 0.1 .* (1 .+ sin.(0.2 .* X)), stn),
]

function compare(solvers, cost)
    for (pname, p) in PROBLEMS
        θ0 = initial_log_θ(p)
        @printf("\n%s  (%d params, start cost %.4g)\n", pname, length(θ0), cost(p, θ0))
        @printf("  %-40s %10s %9s %8s %6s %12s %10s\n", "solver", "median", "memory", "allocs", "iters", "final cost", "‖fit-y‖")
        for (name, solve) in solvers
            try
                θ, iters = solve(p, θ0)  # warm-up / compile
                b = @benchmark $solve($p, $θ0) seconds = 5 evals = 1
                t = median(b)
                @printf("  %-40s %10s %9s %8d %6s %12.5g %10.4g\n", name,
                        BenchmarkTools.prettytime(time(t)), BenchmarkTools.prettymemory(memory(t)),
                        allocs(t), iters === missing ? "-" : string(iters), cost(p, θ), fit_resnorm(p, θ))
            catch e
                @printf("  %-40s FAILED: %s\n", name, first(sprint(showerror, e), 120))
            end
        end
    end
end

println("Julia ", VERSION, ", threads ", Threads.nthreads(), ", maxiter ", MAXITER)
println("\n=== Least-squares solvers (objective: get_lm_objective_func) ===")
compare(LS_SOLVERS, ls_cost)
println("\n=== Quasi-Newton solvers (objective: get_newton_objective_func) ===")
compare(QN_SOLVERS, qn_cost)
