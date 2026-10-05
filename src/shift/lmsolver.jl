# Levenberg-Marquart with a Jacobian that is computed by ForwardDiff only over the
# first np parameters and is constant in the remaining ones. Used for linear
# backgrounds, whose Jacobian columns are fixed by the background basis.
const DiffResults = ForwardDiff.DiffResults
using OptimizationAlgorithms: AbstractGaussNewton, update_direction!, objective!

struct SplitJacobianLM{T, F, FP, D, R, RP, C, JT, BT} <: AbstractGaussNewton{T}
    f::F       # full residual f(rp, θ)
    fp::FP     # f restricted to the first np parameters, the rest read from θ_rest
    d::D       # pre-allocation for parameter direction vector
    r::R       # DiffResult of the full problem
    rp::RP     # DiffResult of the ForwardDiff part
    cfg::C     # Jacobian config of the ForwardDiff part
    JJ::JT     # pre-allocation for Hessian approximation
    J_rest::BT # constant Jacobian columns of the remaining parameters
    θ_rest::Vector{T}
    np::Int
end

function SplitJacobianLM(f, x::AbstractVector{T}, y::AbstractVector, np::Int, J_rest::AbstractMatrix) where T
    θ_rest = x[np+1:end]
    fp = (r, θp) -> f(r, vcat(θp, θ_rest))
    r = DiffResults.JacobianResult(y, x)
    rp = DiffResults.JacobianResult(y, x[1:np])
    cfg = ForwardDiff.JacobianConfig(fp, y, x[1:np])
    J = DiffResults.jacobian(r)
    SplitJacobianLM(f, fp, similar(x), r, rp, cfg, J'J, J_rest, θ_rest, np)
end

function OptimizationAlgorithms.update_jacobian!(LM::SplitJacobianLM, x::AbstractVector)
    y, J = OptimizationAlgorithms.get_value_jacobian(LM)
    np = LM.np
    LM.θ_rest .= @view x[np+1:end]
    ForwardDiff.jacobian!(LM.rp, LM.fp, y, x[1:np], LM.cfg)
    J[:, 1:np] .= DiffResults.jacobian(LM.rp)
    J[:, np+1:end] .= LM.J_rest
    return y, J
end

# Same algorithm as OptimizationAlgorithms.optimize!(::LevenbergMarquart, ...) and
# lm_backtrack!, which only accept the concrete LevenbergMarquart type.
function lm_loop!(LM::AbstractGaussNewton, x::AbstractVector, y::AbstractVector,
                  stn::LevenbergMarquartSettings, λ::Real = 1e-6,
                  verbose::Union{Val{true}, Val{false}} = Val(false))
    oldx = copy(x)
    val, newval = objective!(LM, y, x), Inf
    for i in 1:stn.max_iter
        OptimizationAlgorithms.update_jacobian!(LM, oldx)
        newval, λ = lm_backtrack_generic!(LM, x, oldx, y, λ, stn)
        verbose isa Val{false} || OptimizationAlgorithms.lm_verbose(i, newval, val, λ, y)
        if OptimizationAlgorithms.converged(stn, newval, val, y, i)
            return x, i
        end
        val = newval
        copy!(oldx, x)
    end
    return x, stn.max_iter
end

function lm_backtrack_generic!(LM::AbstractGaussNewton, x::AbstractVector, oldx::AbstractVector,
                               y::AbstractVector, λ::Real, stn::LevenbergMarquartSettings)
    val = 0.
    for j in 1:stn.max_backtrack
        val, dx = update_direction!(LM, λ, stn.min_diagonal)
        @. x = oldx + dx
        newval = objective!(LM, y, x)
        if !(newval ≥ val || isnan(newval) || any(x -> abs(x) > stn.max_step, dx))
            return newval, max(λ / stn.decrease_factor, stn.min_λ)
        end
        λ = min(λ * stn.increase_factor, stn.max_λ)
    end
    copy!(x, oldx) # if backtracking wasn't successful
    return val, λ
end
