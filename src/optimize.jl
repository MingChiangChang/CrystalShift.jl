using OptimizationAlgorithms
"""
    optimize!

    This is the function that each of the node would call on.
    Try to fit a vector of CrystalPhase to the given spectrum
	Treat θ as a vector only.
"""
function optimize!(phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real = .01, mean_θ::AbstractVector = [1., 1.,.2],
                   std_θ::AbstractVector = [3., .01, 1.];
                   maxiter::Int = 32, regularization::Bool = true)
    θ = get_paramters(phases)
	#TODO Modify simple three priors to fit num of paramter in phases
    optimize!(θ, phases, x, y, std_noise, mean_θ, std_θ,
	          maxiter = maxiter, regularization = regularization)
    for (i, cp) in enumerate(phases)
        phases[i] = CrystalPhase(cp, θ)
		deleteat!(θ, collect(1:get_param_nums(phase[i])))
	end
	return phases
end

# Single phase situation. Put phase into [phase].
function optimize!(phase::CrystalPhase, x::AbstractVector, y::AbstractVector,
                   std_noise::Real = .01, mean_θ::AbstractVector = [1., 1.,.2],
                   std_θ::AbstractVector = [3., .01, 1.],
                   max_iter::Int = 32, regularization::Bool = true)
    optimize!([phase], x, y, std_noise, mean_θ, std_θ, max)
end

# TODO apply prior for different num of free params ....
function optimize!(θ::AbstractVector, phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean_θ::AbstractVector,
                   std_θ::AbstractVector; maxiter::Int = 32,
                   regularization::Bool = true)
	length(phases) == length(mean_θ) == length(std_θ) || error("phases and prior terms should have the same length")
    function residual!(r::AbstractVector, θ::AbstractVector)
        params = θ # make a copy
        @. r = y
        r -= reconstruct!(phases, params, x)
		r ./= sqrt(2) * std_noise # ???
        return r
	end

	# The lower symmetry phases have better fitting power and thus
	# should be punished more by the prior
    function prior!(p::AbstractVector, θ::AbstractVector)
		μ =  log.(mean_θ) # [0, 0, 0]
		σ² = std_θ.^2 # var_a, var_α, var_σ
		@. p = (θ - μ) / 2σ²
	end

	# Regularized cost function
    function f(rp::AbstractVector, θ::AbstractVector)
    	r = @view rp[1:length(y)] # residual term
    	residual!(r, θ)
    	p = @view rp[length(y)+1:end] # prior term
    	prior!(p, θ)
    	return rp
    end

    # TODO initialize parameters
    initialize_θ!(θ)

    @. θ = log(θ) # tramsform to log space for better conditioning

    if regularization
		r = zeros(eltype(θ), length(y) + length(θ)) # Reason??
        LM = LavernbergMarquart(f, r, θ)
	else
		r = zeros(eltype(θ), size(y))
        LM = LavernbergMarquart(residual!, r, θ)
	end
    stn = LevenbergMarquartSettings(min_resnorm = 1e-2, min_res = 1e-3,
						min_decrease = 1e-8, max_iter = maxiter,
						decrease_factor = 7, increase_factor = 10, max_step = 1.0)
	λ = 1e-6
	OptimizationAlgorithms.optimize!(LM, θ, copy(r), stn, λ, Val(false))

	@. θ = exp(θ) # transform back
end
