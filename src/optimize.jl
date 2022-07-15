function fit_phases(phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean_θ::AbstractVector = [1. , 1.,.2],
                   std_θ::AbstractVector = [1., 1., 1.];
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

function fit_amorphous(W::Wildcard, BG::Background, x::AbstractVector, y::AbstractVector,
					std_noise::Real;
					method::OptimizationMethods,
					objective::String = "LS",
					maxiter::Int = 32,
					regularization::Bool = true,
					verbose::Bool = false, tol::Float64 = DEFAULT_TOL)

    pm = PhaseModel(W, BG)
	opt_stn = OptimizationSettings{Float64}(pm, std_noise, [1., 1., 1.], [1., 1., 1.],
											maxiter, regularization,
											method, objective, verbose, tol)

	y ./= maximum(y) * 2
	opt_pm = optimize!(pm, x, y, opt_stn)
	return opt_pm
end

# TODO: allow change in some parameters
function full_optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector,
	std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
	std_θ::AbstractVector = [1., Inf, 5.];
	method::OptimizationMethods, objective::String = "LS",
	regularization::Bool = true,
	loop_num::Int=8,
	peak_shift_iter::Int = 32,
	mod_peak_num::Int = 32,
	peak_mod_mean::AbstractVector = [1.],
	peak_mod_std::AbstractVector = [.5],
	peak_mod_iter::Int=32,
	verbose::Bool = false, tol::Float64 =DEFAULT_TOL)

	# have_bg = !isnothing(pm.background)
	c = pm
	for i in 1:loop_num
		c = optimize!(c, x, y, std_noise, mean_θ, std_θ;
			method=method, objective=objective, maxiter=peak_shift_iter,
			regularization=regularization, verbose=verbose, tol=tol)

		IMs = get_PeakModCP(c, x, mod_peak_num)

		Mod_IMs = optimize!(IMs, x, y, std_noise, peak_mod_mean, peak_mod_std;
					method=bfgs, objective=objective, maxiter=peak_mod_iter,
					regularization=regularization, verbose=verbose, tol=tol)
		# change_peak_int!.(c.CPs, Mod_IMs[1:end-Int64(have_bg)])
		change_peak_int!.(c.CPs, Mod_IMs)
		# change_c!(c.background, Mod_IMs[end-Int64(have_bg)+1:end])
		c = optimize!(c, x, y, std_noise, mean_θ, std_θ;
			method=method, objective=objective, maxiter=peak_shift_iter,
			regularization=regularization, verbose=verbose, tol=tol)
	end
	return c
end

function full_optimize!(cp::AbstractVector{<:CrystalPhase}, x::AbstractVector, y::AbstractVector,
						std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
						std_θ::AbstractVector = [1., 1., 5.];
						method::OptimizationMethods, objective::String = "LS",
						regularization::Bool = true,
						loop_num::Int=8,
						peak_shift_iter::Int = 32,
						mod_peak_num::Int = 32,
						peak_mod_mean::AbstractVector = [1.],
						peak_mod_std::AbstractVector = [.5],
						peak_mod_iter::Int=32,
						verbose::Bool = false, tol::Float64 =DEFAULT_TOL)
    pm = PhaseModel(cp)
	pm = full_optimize!(pm, x, y, std_noise, mean_θ, std_θ;
						method=method, objective=objective,
						regularization=regularization,
						loop_num=loop_num,
						peak_shift_iter=peak_shift_iter,
						mod_peak_num=mod_peak_num,
						peak_mod_mean=peak_mod_mean,
						peak_mod_std=peak_mod_std,
						peak_mod_iter=peak_mod_iter,
						verbose=verbose, tol=tol)
	pm.CPs
end

function full_optimize!(cp::CrystalPhase, x::AbstractVector, y::AbstractVector,
	std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
	std_θ::AbstractVector = [1., 1., 5.];
	method::OptimizationMethods, objective::String = "LS",
	regularization::Bool = true,
	loop_num::Int=8,
	peak_shift_iter::Int = 32,
	mod_peak_num::Int = 32,
	peak_mod_mean::AbstractVector = [1.],
	peak_mod_std::AbstractVector = [.5],
	peak_mod_iter::Int=32,
	verbose::Bool = false, tol::Float64 =DEFAULT_TOL)

	full_optimize!([cp], x, y, std_noise, mean_θ, std_θ;
			method=method, objective=objective,
			regularization=regularization,
			loop_num=loop_num,
			peak_shift_iter=peak_shift_iter,
			mod_peak_num=mod_peak_num,
			peak_mod_mean=peak_mod_mean,
			peak_mod_std=peak_mod_std,
			peak_mod_iter=peak_mod_iter,
			verbose=verbose, tol=tol)
end
"""
    optimize!

    This is the function that each of the node would call on.
    Try to fit a PhaseModel, which comprise a vector of CrystalPhase
	and an optional BackgroundModel to the given spectrum y.

	Return: a PhaseModel object
"""
function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector,
					std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
					std_θ::AbstractVector = [1., Inf, 5.];
					method::OptimizationMethods, objective::String = "LS",
					maxiter::Int = 32,
					regularization::Bool = true,
					verbose::Bool = false, tol::Float64 =DEFAULT_TOL)
	opt_stn = OptimizationSettings{Float64}(pm, std_noise, mean_θ, std_θ,
							maxiter, regularization,
							method, objective, verbose, tol)

	optimize!(pm, x, y, opt_stn)
end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings)
	θ = get_free_params(pm)
	θ = optimize!(θ, pm, x, y, opt_stn)

	pm = reconstruct!(pm, θ)
	return pm
end

function optimize!(θ::AbstractVector, pm::PhaseModel,
				   x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings)
	if eltype(pm.CPs) <: CrystalPhase
	    θ = initialize_activation!(θ, pm, x, y)
	end
    # TODO: Don't take log of profile parameters
	# @views test
	# or . test
	θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .= @views log.(θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)]) # tramsform to log space for better conditioning
	log_θ = θ
	(any(isnan, log_θ) || any(isinf, log_θ)) && throw("any(isinf, θ) = $(any(isinf, θ)), any(isnan, θ) = $(any(isnan, θ))")

	# TODO use Match.jl, or just use multiple dispatch on method?
	if opt_stn.method == LM
		log_θ = lm_optimize!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == Newton
		log_θ = newton!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == bfgs
		log_θ = BFGS!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == l_bfgs
		log_θ = LBFGS!(log_θ, pm, x, y, opt_stn)
	end

	log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] = @views exp.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
	θ = log_θ
	return θ
end

function optimize_with_uncertainty!(pm::PhaseModel, x::AbstractVector, y::AbstractVector,
		std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
		std_θ::AbstractVector = [1., Inf, 5.];
		method::OptimizationMethods, objective::String = "LS",
		maxiter::Int = 32,
		regularization::Bool = true,
		verbose::Bool = false, tol::Float64 = DEFAULT_TOL)
    objective == "LS" || error("Invalid objective!")
	opt_stn = OptimizationSettings{Float64}(pm, std_noise, mean_θ, std_θ,
				maxiter, regularization,
				method, objective, verbose, tol)

	optimize_with_uncertainty!(pm, x, y, opt_stn)
end

function optimize_with_uncertainty!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings)
	θ = get_free_params(pm)
	θ, uncer = optimize_with_uncertainty!(θ, pm, x, y, opt_stn)

	pm = reconstruct!(pm, θ)
	return pm, uncer
end

function optimize_with_uncertainty!(θ::AbstractVector, pm::PhaseModel,
									x::AbstractVector, y::AbstractVector,
									opt_stn::OptimizationSettings)
	if eltype(pm.CPs) <: CrystalPhase
		θ = initialize_activation!(θ, pm, x, y)
	end
	# TODO: Don't take log of profile parameters
	# @views test
	# or . test
	θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .= @views log.(θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)]) # tramsform to log space for better conditioning
	log_θ = θ
	(any(isnan, log_θ) || any(isinf, log_θ)) && throw("any(isinf, θ) = $(any(isinf, θ)), any(isnan, θ) = $(any(isnan, θ))")

	# TODO use Match.jl, or just use multiple dispatch on method?
	if opt_stn.method == LM
		log_θ = lm_optimize!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == Newton
		log_θ = newton!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == bfgs
		log_θ = BFGS!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == l_bfgs
		log_θ = LBFGS!(log_θ, pm, x, y, opt_stn)
	end

	# Background is linear. Hessian is always 0. Need to remove to prevent a weird inexact error
	phase_params = get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)
	_, new_bg = reconstruct_BG!(exp.(log_θ[phase_params+1:end]), pm.background)
	signal = y .- evaluate!(zero(y), new_bg, x)
	phases = PhaseModel(pm.CPs, pm.wildcard, nothing)
	phase_log_θ = log_θ[1:phase_params]

	if opt_stn.method == LM
		f = get_lm_objective_func(phases, x, signal, opt_stn)
		r = zeros(Real, length(y) + phase_params)
		function res(log_θ)
			sum(abs2, f(r, log_θ))
		end
	else
		res = get_newton_objective_func(pm, x, y, opt_stn)
	end

	H = ForwardDiff.hessian(res, phase_log_θ)
	val = res(phase_log_θ)
	display(val)
	uncer = sqrt.(diag(val / (length(x) - length(phase_log_θ)) * inverse(H)))

	log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .= @views exp.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
	θ = log_θ
	return θ, uncer
end


function initialize_activation!(θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector)
    new_θ = copy(θ) # make a copy
	start = 1
	for phase in pm.CPs
        param_num = get_param_nums(phase)
		p = evaluate(phase, θ[start:start+param_num-1], x)
		new_θ[start + param_num - 2 - get_param_nums(phase.profile)] = dot(p, y) / sum(abs2, p)
        start += param_num
	end
	return new_θ
end

function lm_optimize!(log_θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector, 
	                 opt_stn::OptimizationSettings)
	opt_stn.objective == "LS" || error("LM only work with LS for now")

	f = get_lm_objective_func(pm, x, y, opt_stn)
	if opt_stn.regularization
		r = zeros(eltype(log_θ), length(y) + length(log_θ) )
		LM = LevenbergMarquart(f, log_θ, r)
	else
		r = zeros(eltype(log_θ), size(y))
		LM = LevenbergMarquart(f, log_θ, r)
	end

	stn = LevenbergMarquartSettings(min_resnorm = 1e-2, min_res = 1e-3,
						min_decrease = 1e-6, max_iter = opt_stn.maxiter,
						decrease_factor = 7, increase_factor = 10, max_step = 0.1)

	λ = 1e-6
	OptimizationAlgorithms.optimize!(LM, log_θ, copy(r), stn, λ, Val(opt_stn.verbose))#, false)
	return log_θ
end

function get_lm_objective_func(pm::PhaseModel,
							   x::AbstractVector, y::AbstractVector,
							   opt_stn::OptimizationSettings)
	pr = opt_stn.priors

	function residual!(r::AbstractVector, log_θ::AbstractVector)
		_residual!(pm, log_θ, x, y, r, pr.std_noise)
	end

	function prior!(p::AbstractVector, log_θ::AbstractVector)
		_prior(p, log_θ, pr.mean_θ, pr.std_θ)
	end

	# Regularized cost function
	function f(rp::AbstractVector, log_θ::AbstractVector)
		if (any(isinf, log_θ) || any(isnan, log_θ))
			return Inf
		end
		bg_param_num = get_param_nums(pm.background)
		w_param_num = get_param_nums(pm.wildcard)
		θ_cp = log_θ[1:end - w_param_num-bg_param_num]
		θ_w  = log_θ[end - w_param_num-bg_param_num+1 : end-bg_param_num]
		θ_bg = log_θ[end - bg_param_num + 1 : end]
		r = @view rp[1:length(y)] # residual term
		residual!(r, log_θ)
		p = @view rp[length(y)+1:length(y)+get_param_nums(pm.CPs)] # prior term
		prior!(p, θ_cp)
		wp = @view rp[length(y)+get_param_nums(pm.CPs)+1:length(y)+get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)]
		lm_prior!(wp, pm.wildcard, θ_w)
		bg_p = @view rp[length(y)+get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)+1:end]
		lm_prior!(bg_p, pm.background, θ_bg)
		return rp
	end

	opt_stn.regularization ? (return f) : (return residual!)
end

# TODO: Work with newton
function newton!(log_θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector,
				opt_stn::OptimizationSettings)
	tol, maxiter, verbose = opt_stn.tol, opt_stn.maxiter, opt_stn.verbose

	N = SaddleFreeNewton(get_newton_objective_func(pm, x, y, opt_stn), log_θ)
	N = UnitDirection(N)
	D = DecreasingStep(N, log_θ)
	# IDEA: D = OptimizationAlgorithms.TrustedDirection(D, maxnorm, maxentry)
	S = StoppingCriterion(log_θ, dx = tol, rx = tol,
							maxiter = maxiter, verbose = verbose)
	fixedpoint!(D, log_θ, S)

	return log_θ
end
using OptimizationAlgorithms: UnitDirection
function LBFGS!(log_θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector,
				opt_stn::OptimizationSettings)
	tol, maxiter, verbose = opt_stn.tol, opt_stn.maxiter, opt_stn.verbose

	N = LBFGS(get_newton_objective_func(pm, x, y, opt_stn), log_θ, 10) # default to 10
	N = UnitDirection(N)
	D = DecreasingStep(N, log_θ)
	S = StoppingCriterion(log_θ, dx = tol, rx=tol, maxiter=maxiter, verbose=verbose)
	fixedpoint!(D, log_θ, S)
	return log_θ
end

function BFGS!(log_θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector,
				opt_stn::OptimizationSettings)
	tol, maxiter, verbose = opt_stn.tol, opt_stn.maxiter, opt_stn.verbose

	N = BFGS(get_newton_objective_func(pm, x, y, opt_stn), log_θ)
	N = UnitDirection(N)
	D = DecreasingStep(N, log_θ)
	S = StoppingCriterion(log_θ, dx = tol, rx=tol, maxiter=maxiter, verbose=verbose)
	fixedpoint!(D, log_θ, S)
	return log_θ
end

function get_newton_objective_func(pm::PhaseModel,
									x::AbstractVector, y::AbstractVector,
									opt_stn::OptimizationSettings)
	pr = opt_stn.priors
	mean_log_θ = log.(pr.mean_θ)

	function prior(log_θ::AbstractVector)
		bg_param_num = get_param_nums(pm.background)
		w_param_num = get_param_nums(pm.wildcard)
		θ_cp = log_θ[1:end - w_param_num-bg_param_num]
		θ_w  = log_θ[end - w_param_num-bg_param_num+1 : end-bg_param_num]
		θ_bg = log_θ[end - bg_param_num + 1 : end]
		p = zero(eltype(log_θ))
		@inbounds @simd for i in eachindex(θ_cp)
			p += ((θ_cp[i] - mean_log_θ[i]) / (sqrt(2)*pr.std_θ[i]))^2
		end
		p += _prior(pm.background, θ_bg)
		p += _prior(pm.wildcard, θ_w)
		return p
	end

	# Regularized cost function
	# NOTE on order of inputs in KL divergence:
	# kl(y, r_θ) is more inclusive, i.e. it tries to fit all peaks, even if it can't
	# kl(r_θ, y) is more exclusive, i.e. it tends to fit peaks well that it can explain while ignoring others

	function kl_objective(log_θ::AbstractVector) # TODO: Fix this for wildcard
		log_θ[1:get_param_nums(pm.CPs)] .= @views exp.(log_θ[1:get_param_nums(pm.CPs)])
		if (any(isinf, log_θ) || any(isnan, log_θ))
			log_θ[1:get_param_nums(pm.CPs)] .= @views log.(log_θ[1:get_param_nums(pm.CPs)])
			return Inf
		end
		r_θ = evaluate(pm, log_θ, x) # reconstruction of phases, TODO: pre-allocate result (one for Dual, one for Float)
		r_θ ./= exp(1) # since we are not normalizing the inputs, this rescaling has the effect that kl(α*y, y) has the optimum at α = 1
		log_θ[1:get_param_nums(pm.CPs)] .= @views log.(log_θ[1:get_param_nums(pm.CPs)])
		p_θ = prior(log_θ)
		λ = 0 #TODO: Fix the prior optimization problem and add it to the setting
		kl(r_θ, y) + λ * p_θ
	end

	function ls_residual(log_θ::AbstractVector)
		(any(isinf, log_θ) || any(isnan, log_θ)) && return Inf
		r = zeros(promote_type(eltype(log_θ), eltype(x), eltype(y)), length(x))
		r = _residual!(pm, log_θ, x, y, r, pr.std_noise)
		return sum(abs2, r)
	end

	function ls_prior(log_θ::AbstractVector)
		bg_param_num = get_param_nums(pm.background)
		w_param_num = get_param_nums(pm.wildcard)
		θ_cp = log_θ[1:end - w_param_num-bg_param_num]
		θ_w  = log_θ[end - w_param_num-bg_param_num+1 : end-bg_param_num]
		θ_bg = log_θ[end - bg_param_num + 1 : end]
		p = zero(θ_cp)
		return (sum(abs2, _prior(p, θ_cp, pr.mean_θ, pr.std_θ))
		        + _prior(pm.background, θ_bg)
				+ _prior(pm.wildcard, θ_w) )
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

function optimize!(phases::AbstractVector{<:AbstractPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
                   std_θ::AbstractVector = [1., Inf, 5.];
                   method::OptimizationMethods, objective::String = "LS",
				   maxiter::Int = 32,
				   regularization::Bool = true,
				   verbose::Bool = false, tol::Float64 =DEFAULT_TOL)
	pm = PhaseModel(phases)
	pm = optimize!(pm, x, y, std_noise, mean_θ, std_θ, method=method,
	             objective=objective, maxiter= maxiter,regularization=regularization,
				 verbose=verbose, tol=tol)
    pm.CPs
end


# Single phase situation. Put phase into [phase].
function optimize!(phase::AbstractPhase,
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

############################# Objective Helper ############################
function _prior(p::AbstractVector, log_θ::AbstractVector,
				mean_θ::AbstractVector, std_θ::AbstractVector)
	mean_log_θ = log.(mean_θ)
	# eltype(θ) <: AbstractFloat ?  mean_log_θ : mean_log_θ_dual
	@. p = (log_θ - mean_log_θ) / (sqrt(2)*std_θ)
	return p # IDEA: Is this too small??
end

function _residual!(pm::PhaseModel,
					log_θ::AbstractVector,
					x::AbstractVector, y::AbstractVector,
					r::AbstractVector,
					std_noise::Real)
	log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .= @views exp.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
	(any(isinf, log_θ) || any(isnan, log_θ)) && return Inf
	@. r = y
	evaluate_residual!(pm, log_θ, x, r) # Avoid allocation, put everything in here??
	r ./= sqrt(2) * std_noise # trade-off between prior and
	# actual residual
	log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .=  @views log.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
	return r
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
