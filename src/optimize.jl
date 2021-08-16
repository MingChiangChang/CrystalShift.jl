using OptimizationAlgorithms
using OptimizationAlgorithms: LevenbergMarquart, LevenbergMarquartSettings
using OptimizationAlgorithms: update_jacobian!
using LinearAlgebra
using ForwardDiff: Dual
# IDEA might want to define a prior object

"""
    optimize!

    This is the function that each of the node would call on.
    Try to fit a vector of CrystalPhase to the given spectrum y.
	Treat θ as a vector only.
"""
function optimize!(phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real = .01, mean_θ::AbstractVector = [1., 1.,.2],
                   std_θ::AbstractVector = [3., .01, 1.];
                   maxiter::Int = 32, regularization::Bool = true)
    θ = get_parameters(phases)
    println(θ)
	if length(mean_θ) == 3
	    mean_θ, std_θ = extend_priors(mean_θ, std_θ, phases)
	end
    println(mean_θ, std_θ)
	(sum([get_param_nums(phase) for phase in phases]) == length(mean_θ) == length(std_θ) ||
	 error("number of parameter must match number of terms in the prior"))

    θ = optimize!(θ, phases, x, y, std_noise, mean_θ, std_θ,
	          maxiter = maxiter, regularization = regularization)

    for (i, cp) in enumerate(phases)
        phases[i] = CrystalPhase(cp, θ)
		deleteat!(θ, collect(1:get_param_nums(phases[i])))
	end
	return phases
end

function extend_priors(mean_θ::AbstractVector, std_θ::AbstractVector,
	                    phases::AbstractVector{<:CrystalPhase})
	totl_params = sum([get_param_nums(phase) for phase in phases])
	full_mean_θ = zeros(totl_params)
	full_std_θ = zeros(totl_params)
	start = 1
	for phase in phases
		n = phase.cl.free_param
		println(n)
		full_mean_θ[start:start+n-1] = repeat(mean_θ[1, :], n)
		full_std_θ[start:start+n-1] = repeat(std_θ[1, :], n)
		full_mean_θ[start + n: start + n + 1] = mean_θ[2:3]
		full_std_θ[start + n:start + n + 1] = std_θ[2:3]
		start += (n+2)
    end
	return full_mean_θ, full_std_θ
end

# Single phase situation. Put phase into [phase].
function optimize!(phase::CrystalPhase, x::AbstractVector, y::AbstractVector,
                   std_noise::Real = .01, mean_θ::AbstractVector = [1., 1.,.2],
                   std_θ::AbstractVector = [.1, .5, 1.];
                   maxiter::Int = 32, regularization::Bool = true)
	optimize!([phase], x, y, std_noise, mean_θ, std_θ,
	          maxiter=maxiter, regularization=regularization)
end

function initialize_activation!(θ::AbstractVector, phases::AbstractVector,
	                 x::AbstractVector, y::AbstractVector)
	# Use inner product to set initial activation
	new_θ = copy(θ) # make a copy
	start = 1
	for phase in phases
        param_num = get_param_nums(phase)
		p = reconstruct!(phase, θ, x)
		new_θ[start + param_num - 2] =  dot(p, y) / sum(abs2, p)
        start += param_num
	end
	return new_θ
end

# TODO apply prior for different num of free params
function optimize!(θ::AbstractVector, phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean_θ::AbstractVector,
                   std_θ::AbstractVector; maxiter::Int = 32,
                   regularization::Bool = true)
    function residual!(r::AbstractVector, θ::AbstractVector)
        params = exp.(θ) # make a copy
        @. r = y
		r .-= reconstruct!(phases, params, x)
		#r ./= sqrt(2) * std_noise # ???
		# if r isa AbstractVector{<:Dual}
		# 	plot_r = [r[i].value for i in eachindex(r)]
		# 	plt = plot(x, plot_r, title="Residual")
		# else
		#     plt = plot(x, r, title="Residual")
		# end
		# display(plt)
		#savefig("test_$(sum(r)).png")
        return r
	end

	# The lower symmetry phases have better fitting power and thus
	# should be punished more by the prior
    function prior!(p::AbstractVector, θ::AbstractVector)
		μ = log.(mean_θ) # [0, 0, 0]
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
    θ = initialize_activation!(θ, phases, x, y)

    @. θ = log(θ) # tramsform to log space for better conditioning
    (any(isnan, θ) || any(isinf, θ)) && throw("any(isinf, θ) = $(any(isinf, θ)), any(isnan, θ) = $(any(isnan, θ))")

    if regularization
		r = zeros(eltype(θ), length(y) + length(θ)) # Reason??
        LM = LevenbergMarquart(f, θ, r)
	else
		r = zeros(eltype(θ), size(y))
        LM = LevenbergMarquart(residual!, θ, r)
	end
    stn = LevenbergMarquartSettings(min_resnorm = 1e-2, min_res = 1e-3,
						min_decrease = 1e-8, max_iter = maxiter,
						decrease_factor = 7, increase_factor = 10, max_step = 1.0)
	λ = 1e-6
	#t, J = update_jacobian!(LM, θ)

	OptimizationAlgorithms.optimize!(LM, θ, copy(r), stn, λ, Val(false))
	@. θ = exp(θ) # transform back
	return θ
end
