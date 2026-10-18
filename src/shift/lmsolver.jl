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
