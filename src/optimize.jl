using OptimizationAlgorithms
"""
    optimize!

    This is the function that each of the node would call on.
    Try to fit a vector of CrystalPhase to the given spectrum
"""
function optimize!(phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Reak = .01, mean_θ::AbstractVector = [1., 1.,.2],
                   std_θ::AbstractVector = [3., .01, 1.];
                   max_iter::Int = 32, regularization::Bool = true)
    return
end

# Single phase situation. Put phase into [phase].
function optimize!(phase::CrystalPhase, x::AbstractVector, y::AbstractVector,
                   std_noise::Reak = .01, mean_θ::AbstractVector = [1., 1.,.2],
                   std_θ::AbstractVector = [3., .01, 1.],
                   max_iter::Int = 32, regularization::Bool = true)
    optimize!([phase], x, y, std_noise, mean_θ, std_θ, max)
end

# TODO Line 62 in PhaseMapping/src/optimization.jl

function optimize!(θ::AbstracVector, phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean::AbstractVector,
                   std_θ::AbstractVector; maxiter::Int = 32,
                   regularization::Bool = true)
    function residual!(r::AbstractVector, θ::AbstractVector)
        params = θ # make a copy
        @. r = y
        r -= reconstruct!(phases, θ, x)
        return r
    # l2 prior
    function prior!(p::AbstractVector, θ::AbstractVector) end
    function f(rp::AbstractVector, θ::AbstractVector)
    	r = @view rp[1:length(y)] # residual term
    	residual!(r, θ)
    	p = @view rp[length(y)+1:end] # prior term
    	prior!(p, θ)
    	return rp
    end

    if regularization
        LM = LavernbergMarquart(f, r, θ)
	else
        LM = LavernbergMarquart(residual!, r, θ)
	end
    stn = LevenbergMarquartSettings(min_resnorm = 1e-2, min_res = 1e-3,
						min_decrease = 1e-8, max_iter = maxiter,
						decrease_factor = 7, increase_factor = 10, max_step = 1.0)
						λ = 1e-6
	OptimizationAlgorithms.optimize!(LM, θ, copy(r), stn, λ, Val(false))

	return θ
end
