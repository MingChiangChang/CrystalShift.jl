# Levenberg-Marquart with a Jacobian assembled from independent column blocks. The LM
# residual (see get_lm_objective_func) is w .* (y - sum of component models) in the data
# rows, plus one prior row per parameter (row n + i for parameter i). Each phase and the
# wildcard only depend on their own parameters, so ForwardDiff runs over each of them on
# its own and its Dual numbers carry only that component's partials instead of all of
# them. A linear background (is_linear trait) contributes constant Jacobian columns from
# its basis. The residual value is assembled from the same block evaluations.
# CrystalPhase blocks with a supported profile use a hand-written Jacobian instead
# (PhaseJacobianBlock, analyticjacobian.jl), falling back to their ForwardDiff block.
const DiffResults = ForwardDiff.DiffResults
using OptimizationAlgorithms: AbstractGaussNewton, update_direction!, objective!

# A block has cols, out and J, and block_jacobian!(b, θ_block) fills out and J
abstract type JacobianBlock end

# g!(out, θ_block) writes the block's data rows (n) followed by its prior rows, if any
struct ADBlock{G, C} <: JacobianBlock
    g!::G
    cols::UnitRange{Int}
    out::Vector{Float64}
    J::Matrix{Float64}
    cfg::C
end

function ADBlock(g!, cols::UnitRange{Int}, rows::Int)
    out = zeros(rows)
    ADBlock(g!, cols, out, zeros(rows, length(cols)), ForwardDiff.JacobianConfig(g!, out, zeros(length(cols))))
end

# value!(out, θ_block) writes the block's data rows followed by its prior rows; J is constant
struct ConstantBlock{V}
    value!::V
    cols::UnitRange{Int}
    out::Vector{Float64}
    J::Matrix{Float64}
end

struct BlockJacobianLM{T, F, D, R, JT, B, C} <: AbstractGaussNewton{T}
    f::F            # full residual f(rp, θ), used for objective values
    d::D            # pre-allocation for parameter direction vector
    r::R            # DiffResult holding the assembled residual and Jacobian
    JJ::JT          # pre-allocation for Hessian approximation
    blocks::B       # JacobianBlocks, one per phase and one for the wildcard
    constant::C     # ConstantBlock of a linear background, or nothing
    base::Vector{T} # w .* y in the data rows, zero in the prior rows
    n::Int          # number of data rows
end

function BlockJacobianLM(f, x::AbstractVector{T}, y::AbstractVector, blocks, constant,
                         base::AbstractVector, n::Int) where T
    r = DiffResults.JacobianResult(y, x)
    J = DiffResults.jacobian(r)
    BlockJacobianLM(f, similar(x), r, J'J, blocks, constant, Vector{T}(base), n)
end

function add_block!(y::AbstractVector, J::AbstractMatrix, b, n::Int)
    @views y[1:n] .+= b.out[1:n]
    @views J[1:n, b.cols] .= b.J[1:n, :]
    if length(b.out) > n # prior rows
        @views y[n .+ b.cols] .+= b.out[n+1:end]
        @views J[n .+ b.cols, b.cols] .= b.J[n+1:end, :]
    end
end

function OptimizationAlgorithms.update_jacobian!(LM::BlockJacobianLM, x::AbstractVector)
    y, J = OptimizationAlgorithms.get_value_jacobian(LM)
    y .= LM.base
    J .= 0
    for b in LM.blocks
        block_jacobian!(b, x[b.cols])
        add_block!(y, J, b, LM.n)
    end
    if !isnothing(LM.constant)
        b = LM.constant
        b.value!(b.out, x[b.cols])
        add_block!(y, J, b, LM.n)
    end
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

# Powell's dogleg with a scaled trust region on the same BlockJacobianLM, so it uses the block
# (and hand-written) Jacobians. Same algorithm, constants and box handling as LeastSquaresOptim's
# Dogleg (v1, optimizer/dogleg.jl), except:
#  - the Gauss-Newton step solves the normal equations J'J δ = J'r (Cholesky of the p×p J'J,
#    falling back to QR if it is not positive definite) instead of a QR of the (n+p)×p Jacobian;
#  - a rejected Gauss-Newton step that would still fit in the halved trust region is not retried:
#    with x and J unchanged the retry is identical and is rejected again, so Δ is halved past
#    the step at once. Each skipped halving still counts as an iteration, so steps, Δ and the
#    result are exactly those of retrying, without the residual evaluations;
#  - it also stops when the step's predicted decrease of ‖r‖² is below pred_tol·‖r‖²
#    (DOGLEG_PRED_TOL = 1e-6, about LM's min_decrease; 0 turns it off). Near the minimum the
#    actual decrease is round-off, so steps are rejected until the 1e-8 x tolerance stops it.
# Minimizes ‖r(x)‖² with r = LM.f; lower/upper are box bounds on x (empty for none).
# Returns (x, iterations, converged); x is left at the last accepted point.
const DOGLEG_MIN_STEP_QUALITY = 1e-3
const DOGLEG_MIN_Δ = 1e-16
const DOGLEG_MIN_DIAGONAL, DOGLEG_MAX_DIAGONAL = 1e-6, 1e32
const DOGLEG_PRED_TOL = 1e-6

function dogleg_loop!(LM::BlockJacobianLM, x::AbstractVector, r::AbstractVector,
                      lower::AbstractVector, upper::AbstractVector;
                      iterations::Integer = 1000, x_tol::Real = 1e-8, f_tol::Real = 1e-8,
                      g_tol::Real = 1e-8, Δ::Real = 1.0, pred_tol::Real = DOGLEG_PRED_TOL)
    residual!(out, θ) = (v = LM.f(out, θ); v isa Number && fill!(out, Inf); out) # scalar Inf: invalid θ
    p = length(x)
    δgn, δgr, δx, dtd = zeros(p), zeros(p), zeros(p), zeros(p)
    fcur, ftrial, fpredict = copy(r), similar(r), similar(r)
    JJ = LM.JJ
    wnorm(v) = sqrt(sum(i -> dtd[i] * v[i]^2, eachindex(v)))
    residual!(fcur, x)
    ssr = sum(abs2, fcur)
    J = OptimizationAlgorithms.get_value_jacobian(LM)[2]
    reuse, converged, iter = false, false, 0
    wnorm_δgn = wnorm_δgr = α = 0.0
    maxabs_gr = Inf
    while !converged && iter < iterations
        iter += 1
        all(isfinite, x) || break
        if !reuse
            OptimizationAlgorithms.update_jacobian!(LM, x)
            for j in 1:p
                dtd[j] = clamp(sum(abs2, view(J, :, j)), DOGLEG_MIN_DIAGONAL, DOGLEG_MAX_DIAGONAL)
            end
            iter == 1 && (w = wnorm(x); w > 0 && (Δ *= w))
            mul!(δgr, J', fcur) # gradient
            maxabs_gr = maxabs_projected_gradient(δgr, x, lower, upper)
            δgr ./= dtd         # steepest descent in the dtd-scaled metric
            wnorm_δgr = wnorm(δgr)
            mul!(fpredict, J, δgr)
            α = wnorm_δgr^2 / sum(abs2, fpredict) # Cauchy point α δgr
            mul!(JJ, J', J)     # Gauss-Newton step: J'J δgn = J'r
            mul!(δgn, J', fcur)
            F = cholesky!(Symmetric(JJ), check = false)
            if issuccess(F)
                ldiv!(F, δgn)
            else
                δgn .= qr(J, ColumnNorm()) \ fcur
            end
            wnorm_δgn = wnorm(δgn)
        end
        if wnorm_δgn <= Δ          # Gauss-Newton step inside the trust region
            copyto!(δx, δgn)
            wnorm_δx = wnorm_δgn
        elseif wnorm_δgr * α >= Δ  # Cauchy point outside: scaled steepest descent
            δx .= δgr .* (Δ / wnorm_δgr)
            wnorm_δx = Δ
        else                       # dogleg between the Cauchy point and Gauss-Newton step
            b_dot_a = α * sum(i -> dtd[i] * δgr[i] * δgn[i], 1:p)
            a_sq = (α * wnorm_δgr)^2
            b_minus_a_sq = a_sq - 2b_dot_a + wnorm_δgn^2
            c = b_dot_a - a_sq
            d = sqrt(c^2 + b_minus_a_sq * (Δ^2 - a_sq))
            β = c <= 0 ? (d - c) / b_minus_a_sq : (Δ^2 - a_sq) / (d + c)
            @. δx = β * δgn + α * (1 - β) * δgr
            wnorm_δx = wnorm(δx)
        end
        isempty(lower) || (@. δx = min(δx, x - lower)) # box constraints
        isempty(upper) || (@. δx = max(δx, x - upper))
        x .-= δx
        residual!(ftrial, x)
        trial_ssr = sum(abs2, ftrial)
        mul!(fpredict, J, δx)
        fpredict .-= fcur
        predicted_reduction = abs(ssr - sum(abs2, fpredict))
        ρ = predicted_reduction > 0 ? (ssr - trial_ssr) / predicted_reduction : zero(ssr)
        accepted = ρ >= DOGLEG_MIN_STEP_QUALITY
        converged = (accepted && abs(trial_ssr - ssr) <= f_tol * (abs(ssr) + f_tol)) ||
                    maximum(abs, δx) <= x_tol || maxabs_gr <= g_tol ||
                    predicted_reduction <= pred_tol * ssr
        if accepted
            reuse = false
            copyto!(fcur, ftrial)
            ssr = trial_ssr
        else
            reuse = true
            x .+= δx
        end
        if ρ < 0.25
            Δ = max(DOGLEG_MIN_Δ, Δ / 2)
            if !accepted && wnorm_δx == wnorm_δgn # rejected Gauss-Newton step: skip identical retries
                while !converged && wnorm_δgn <= Δ && Δ > DOGLEG_MIN_Δ && iter < iterations
                    Δ = max(DOGLEG_MIN_Δ, Δ / 2)
                    iter += 1
                end
            end
        elseif ρ > 0.75
            Δ = max(Δ, 3wnorm_δx)
        end
    end
    copyto!(r, fcur)
    return x, iter, converged
end

# ∞-norm of the gradient with components that push into an active bound set to zero
function maxabs_projected_gradient(g::AbstractVector, x::AbstractVector, lower::AbstractVector, upper::AbstractVector)
    m = 0.0
    for i in eachindex(g)
        gi = g[i]
        if !isempty(lower) && x[i] <= lower[i] && gi > 0
            gi = zero(gi)
        elseif !isempty(upper) && x[i] >= upper[i] && gi < 0
            gi = zero(gi)
        end
        m = max(m, abs(gi))
    end
    m
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
