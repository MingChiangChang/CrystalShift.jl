function fit_phases(phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean_θ::AbstractVector = [1. ,.2],
                   std_θ::AbstractVector = [1., 1.];
                   maxiter::Int = 32, regularization::Bool = true)
	optimized_phases = Vector{CrystalPhase}(undef, size(phases))
    @threads for i in eachindex(phases)
        optimized_phases[i] = optimize!(phases[i], x, y, std_noise, mean_θ, std_θ,
	                        maxiter=maxiter, regularization=regularization)[1]
    end
	return optimized_phases[get_min_index(optimized_phases, x, y)]
end

function get_min_index(optimized_phases::AbstractVector{<:CrystalPhase},
	                   x::AbstractVector, y::AbstractVector)
   argmin([norm(p.(x)-y) for p in optimized_phases])
end

"""
    optimize!

    This is the function that each of the node would call on.
    Try to fit a vector of CrystalPhase to the given spectrum y.
	Treat θ as a vector only.
"""
function optimize!(phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
                   std_θ::AbstractVector = [1., Inf, 5.];
                   method::OptimizationMethods, objective::String = "LS",
				   maxiter::Int = 32,
				   regularization::Bool = true,
				   verbose::Bool = false, tol::Float64 =DEFAULT_TOL)  
				   
    opt_stn = OptimizationSettings{Float64}(phases, std_noise, mean_θ, std_θ, 
										    maxiter, regularization, 
										    method, objective, verbose, tol)

    θ = optimize!(phases, x, y, opt_stn)

	start = 1
	for (i, cp) in enumerate(phases)
		phases[i] = CrystalPhase(cp, θ[start:start + get_param_nums(phases[i])-1])
		start += get_param_nums(phases[i])
	end

	return phases
end

function optimize!(phases::AbstractVector{<:CrystalPhase},
	               x::AbstractVector, y::AbstractVector,
                   opt_stn::OptimizationSettings)
	θ = get_free_params(phases)
	optimize!(θ, phases, x, y, opt_stn)
end

# Single phase situation. Put phase into [phase].
function optimize!(phase::CrystalPhase,
					x::AbstractVector, y::AbstractVector,
					std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
					std_θ::AbstractVector = [1., Inf, 5.];
					method::OptimizationMethods, objective::String = "LS",
					maxiter::Int = 32,
					regularization::Bool = true,
					verbose::Bool = false, tol::Float64 =DEFAULT_TOL)

    optimize!([phase], x, y, std_noise, mean_θ, std_θ,
               method=method, objective= objective, maxiter=maxiter, regularization=regularization,
			   verbose=verbose, tol=tol)
end

function length_check(phases::AbstractVector, mean_θ::AbstractVector, std_θ::AbstractVector)
    (sum([phase.cl.free_param + 2 for phase in phases]) == length(mean_θ) == length(std_θ) )
end

function initialize_activation!(θ::AbstractVector, phases::AbstractVector,
	                 x::AbstractVector, y::AbstractVector, min_activation::Real = 1e-4)
	# Use inner product to set initial activation
	new_θ = copy(θ) # make a copy
	start = 1
	for phase in phases
        param_num = get_param_nums(phase)
		p = evaluate(phase, θ, x)
		new_θ[start + param_num - 2 - get_param_nums(phase.profile)] = dot(p, y) / sum(abs2, p)
        start += param_num
	end
	return new_θ
end

function optimize!(θ::AbstractVector, phases::AbstractVector{<:CrystalPhase},
	               x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings)
    θ = initialize_activation!(θ, phases, x, y)

    @. θ = log(θ) # tramsform to log space for better conditioning
	log_θ = θ
    (any(isnan, log_θ) || any(isinf, log_θ)) && throw("any(isinf, θ) = $(any(isinf, θ)), any(isnan, θ) = $(any(isnan, θ))")

	# TODO use Match.jl
    if opt_stn.method == LM
        log_θ = lm_optimize!(log_θ, phases, x, y, opt_stn)
	elseif opt_stn.method == Newton
        log_θ = newton!(log_θ, phases, x, y, opt_stn)
	end

    @. θ = exp(log_θ)
	return θ
end

function lm_optimize!(log_θ::AbstractVector, phases::AbstractVector{<:CrystalPhase},
	                  x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings)
    opt_stn.objective == "LS" || error("LM only work with LS for now")
	
    f = get_lm_objective_func(phases, x, y, opt_stn)
	if opt_stn.regularization
		r = zeros(eltype(log_θ), length(y) + length(log_θ) ) # Reason?? - size(phases, 1)
        LM = LevenbergMarquart(f, log_θ, r)
	else
		r = zeros(eltype(log_θ), size(y))
        LM = LevenbergMarquart(f, log_θ, r)
	end

	stn = LevenbergMarquartSettings(min_resnorm = 1e-2, min_res = 1e-3,
						min_decrease = 1e-8, max_iter = opt_stn.maxiter,
						decrease_factor = 7, increase_factor = 10, max_step = 0.1)

	λ = 1e-6
	OptimizationAlgorithms.optimize!(LM, log_θ, copy(r), stn, λ, Val(opt_stn.verbose))
	return log_θ
end

function get_lm_objective_func(phases::AbstractVector{<:CrystalPhase},
							x::AbstractVector, y::AbstractVector,
							opt_stn::OptimizationSettings) 
    pr = opt_stn.priors
	function residual!(r::AbstractVector, log_θ::AbstractVector)
		_residual!(phases, log_θ, x, y, r, pr.std_noise)
	end

	function prior!(p::AbstractVector, log_θ::AbstractVector)
		_prior(p, log_θ, pr.mean_θ, pr.std_θ)
	end

	# Regularized cost function
	function f(rp::AbstractVector, log_θ::AbstractVector)
		r = @view rp[1:length(y)] # residual term
		residual!(r, log_θ)
		p = @view rp[length(y)+1:end] # prior term
		prior!(p, log_θ)
		return rp
	end

	opt_stn.regularization ? (return f) : (return residual!)
end 

# optimization based on SaddleFreeNewton method
function newton!(log_θ::AbstractVector, phases::AbstractVector{<:CrystalPhase},
                 x::AbstractVector, y::AbstractVector,
				 opt_stn::OptimizationSettings)
	
	tol, maxiter, verbose = opt_stn.tol, opt_stn.maxiter, opt_stn.verbose

	N = SaddleFreeNewton(get_newton_objective_func(phases, x, y, opt_stn), log_θ)
	D = DecreasingStep(N, log_θ)
	# IDEA: D = OptimizationAlgorithms.TrustedDirection(D, maxnorm, maxentry)	
	S = StoppingCriterion(log_θ, dx = tol, rx = tol,
						  maxiter = maxiter, verbose = verbose)
	fixedpoint!(D, log_θ, S)

	return log_θ
end

function get_newton_objective_func(phases::AbstractVector{<:CrystalPhase},
									x::AbstractVector, y::AbstractVector,
									opt_stn::OptimizationSettings)
	pr = opt_stn.priors
	mean_log_θ = log.(pr.mean_θ)

    function prior(log_θ::AbstractVector)
		p = zero(eltype(log_θ))
		@inbounds @simd for i in eachindex(log_θ)
			p += ((log_θ[i] - mean_log_θ[i]) / (sqrt(2)*pr.std_θ[i]))^2
		end
		return p
	end

	# Regularized cost function
	# NOTE on order of inputs in KL divergence:
	# kl(y, r_θ) is more inclusive, i.e. it tries to fit all peaks, even if it can't
	# kl(r_θ, y) is more exclusive, i.e. it tends to fit peaks well that it can explain while ignoring others
    function kl_objective(log_θ::AbstractVector)
		θ = exp.(log_θ)
		r_θ = evaluate(phases, θ, x) # reconstruction of phases, IDEA: pre-allocate result (one for Dual, one for Float)
		r_θ ./= exp(1) # since we are not normalizing the inputs, this rescaling has the effect that kl(α*y, y) has the optimum at α = 1
		p_θ = prior(log_θ)
		λ = 0 #TODO: Fix the prior optimization problem and add it to the setting
		kl(r_θ, y) + λ * p_θ
    end

	function ls_residual(log_θ::AbstractVector)
		r = zeros(promote_type(eltype(log_θ), eltype(x), eltype(y)), length(x))
		r = _residual!(phases, log_θ, x, y, r, pr.std_noise)
		return sum(abs2, r)
	end
	
	function ls_prior(log_θ::AbstractVector)
		p = zero(log_θ)
		sum(abs2, _prior(p, log_θ, pr.mean_θ, pr.std_θ))
	end

	function ls_objective(log_θ::AbstractVector)
		ls_residual(log_θ) + ls_prior(log_θ)
	end

	if opt_stn.objective == "KL"
		return kl_objective
	elseif opt_stn.objective == "LS"
		return ls_objective
	end
end

############################# Objective Helper ############################
function _residual!(phases::AbstractVector{<:CrystalPhase},
					log_θ::AbstractVector,
					x::AbstractVector, y::AbstractVector,
					r::AbstractVector,
					std_noise::Real)
	params = exp.(log_θ) # make a copy
	@. r = y
	evaluate_residual!(phases, params, x, r) # Avoid allocation, put everything in here??
	r ./= sqrt(2) * std_noise # trade-off between prior and
					# actual residual
	return r
end

function _prior(p::AbstractVector, log_θ::AbstractVector,
				mean_θ::AbstractVector, std_θ::AbstractVector)
	mean_log_θ = log.(mean_θ)
	@. p = (log_θ - mean_log_θ) / (sqrt(2)*std_θ)
	return p
end

########################### parameter helpers ##################################
function check_objective(objective::String)
	objective in ALLOWED_OBJECTIVE || error("objective $(objective) not a allowed objective string")
end

# NOTE: candidates for removal?
function remove_act_from_θ(θ::AbstractVector,
	                      phases::AbstractVector{<:Crystal})
	θ_c = copy(θ)
    cursor = 0
	for phase in phases
	    deleteat!(θ_c, cursor + phase.free_param + 1)
        cursor += phase.free_param + 1
	end
	θ_c
end

function remove_act_from_θ(θ::AbstractVector,
	                      phases::AbstractVector{<:CrystalPhase})
	remove_act_from_θ(θ, collect_crystals(phases))
end
