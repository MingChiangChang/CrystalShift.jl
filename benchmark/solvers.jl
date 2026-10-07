# Optimization backends wrapped behind one interface, all solving the exact problem that
# CrystalShift's `simple_optimize!` builds (same starting point, log-space parameters,
# objective closure and maxiter). Included by optimizers.jl and scenarios.jl.
using CrystalShift
using CrystalShift: PhaseModel, OptimizationSettings, CrystalPhase, LM, bfgs, l_bfgs,
                    initialize_activation!, get_lm_objective_func, get_newton_objective_func,
                    get_free_params, get_param_nums, reconstruct!, evaluate!, BFGS!, LBFGS!
import OptimizationAlgorithms
using OptimizationAlgorithms: LevenbergMarquart, LevenbergMarquartSettings
import NonlinearSolve as NLS
import LeastSquaresOptim as LSO
import Optim
using ADTypes: AutoForwardDiff
using LinearAlgebra

struct Problem
    pm::PhaseModel
    x::Vector{Float64}
    y::Vector{Float64}
    stn::OptimizationSettings{Float64}
end

nlog(p::Problem) = get_param_nums(p.pm.CPs) + get_param_nums(p.pm.wildcard)

# Starting point exactly as in simple_optimize!
function initial_log_θ(p::Problem)
    θ = get_free_params(p.pm)
    if p.pm.CPs isa AbstractVector{<:CrystalPhase}
        θ = initialize_activation!(θ, p.pm, p.x, p.y)
    end
    n = nlog(p)
    θ[1:n] .= log.(θ[1:n])
    θ
end

# Optimized PhaseModel from a log-space solution, as simple_optimize! returns it
function to_model(p::Problem, log_θ)
    θ = copy(log_θ)
    n = nlog(p)
    θ[1:n] .= exp.(θ[1:n])
    reconstruct!(p.pm, θ)
end

fit_resnorm(p::Problem, log_θ) = norm(evaluate!(zero(p.x), to_model(p, log_θ), p.x) .- p.y)

# CrystalShift's objectives return a scalar Inf on non-finite input; in-place solvers need a vector.
function lm_residual(p::Problem)
    f = get_lm_objective_func(p.pm, p.x, p.y, zero(p.y), p.stn)
    function f!(r, θ)
        out = f(r, θ)
        out isa Number && fill!(r, Inf)
        r
    end
end
nresid(p::Problem) = length(p.y) + get_param_nums(p.pm)

ls_cost(p::Problem, θ) = sum(abs2, lm_residual(p)(zeros(nresid(p)), θ))
qn_cost(p::Problem, θ) = get_newton_objective_func(p.pm, p.x, p.y, p.stn)(θ)

# --- solvers: (problem, log_θ0) -> (log_θ, iterations) ---------------------------------------

# The package's own LM and Dogleg (method = LM / dogleg in optimize!). lm_optimize! uses
# SplitJacobianLM, with the exact background Jacobian, when the background is linear.
# Neither reports its iteration count.
current_lm(p, θ) = (CrystalShift.lm_optimize!(copy(θ), p.pm, p.x, p.y, zero(p.y), p.stn), missing)
package_dogleg(p, θ) = (CrystalShift.dogleg_optimize!(copy(θ), p.pm, p.x, p.y, zero(p.y), p.stn), missing)
current_bfgs(p, θ) = (BFGS!(copy(θ), p.pm, p.x, p.y, p.stn), missing)
current_lbfgs(p, θ) = (LBFGS!(copy(θ), p.pm, p.x, p.y, p.stn), missing)

function nls(alg)
    function (p, θ)
        f! = lm_residual(p)
        fn = NLS.NonlinearFunction{true}((r, u, _) -> f!(r, u); resid_prototype = zeros(nresid(p)))
        sol = NLS.solve(NLS.NonlinearLeastSquaresProblem(fn, copy(θ)), alg; maxiters = p.stn.maxiter)
        sol.u, sol.stats.nsteps
    end
end

function lso(alg)
    function (p, θ)
        prob = LSO.LeastSquaresProblem(x = copy(θ), f! = lm_residual(p), output_length = nresid(p), autodiff = :forward)
        res = LSO.optimize!(prob, alg; iterations = p.stn.maxiter)
        res.minimizer, res.iterations
    end
end

function optim(alg)
    function (p, θ)
        cost = get_newton_objective_func(p.pm, p.x, p.y, p.stn)
        res = Optim.optimize(cost, copy(θ), alg, Optim.Options(iterations = p.stn.maxiter); autodiff = AutoForwardDiff())
        Optim.minimizer(res), Optim.iterations(res)
    end
end

# Least-squares solvers minimize ‖lm residual‖²; quasi-Newton solvers minimize the scalar
# newton objective. Each family is scored on its own objective.
const LS_SOLVERS = [
    "current: CrystalShift LM"             => current_lm,
    "CrystalShift dogleg (method = dogleg)" => package_dogleg,
    "NonlinearSolve LevenbergMarquardt"    => nls(NLS.LevenbergMarquardt(; autodiff = AutoForwardDiff())),
    "NonlinearSolve TrustRegion"           => nls(NLS.TrustRegion(; autodiff = AutoForwardDiff())),
    "NonlinearSolve GaussNewton"           => nls(NLS.GaussNewton(; autodiff = AutoForwardDiff())),
    "LeastSquaresOptim LevenbergMarquardt" => lso(LSO.LevenbergMarquardt()),
    "LeastSquaresOptim Dogleg"             => lso(LSO.Dogleg()),
]
const QN_SOLVERS = [
    "current: OptimizationAlgorithms BFGS"  => current_bfgs,
    "current: OptimizationAlgorithms LBFGS" => current_lbfgs,
    "Optim BFGS"  => optim(Optim.BFGS()),
    "Optim LBFGS" => optim(Optim.LBFGS()),
]
