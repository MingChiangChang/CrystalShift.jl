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
					optimize_mode::OptimizationMode=Simple,
					maxiter::Int = 32,
					regularization::Bool = true,
					em_loop_num::Integer = 8,
					λ::Float64 = 1.,
					verbose::Bool = false, tol::Float64 = DEFAULT_TOL)

    pm = PhaseModel(W, BG)
	opt_stn = OptimizationSettings{Float64}(std_noise, [1., 1., 1.], [1., 1., 1.],
											maxiter, regularization,
											method, objective, optimize_mode, em_loop_num, λ, verbose, tol)

	y ./= maximum(y) * 2
	opt_pm = optimize!(pm, x, y, opt_stn)
	return opt_pm
end

# TODO: allow change in some parameters
function full_optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector,
						std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
						std_θ::AbstractVector = [1., Inf, 5.];
						method::OptimizationMethods = LeastSquares(),
						optimize_mode::OptimizationMode=Simple,
						objective::String = "LS",
						regularization::Bool = true,
						loop_num::Int=8,
						peak_shift_iter::Int = 32,
						mod_peak_num::Int = 32,
						peak_mod_mean::AbstractVector = [1.],
						peak_mod_std::AbstractVector = [.5],
						peak_mod_iter::Int=32,
						λ::Float64=1.,
						verbose::Bool = false, tol::Float64 =DEFAULT_TOL)

	# have_bg = !isnothing(pm.background)
	c = pm
	for i in 1:loop_num
		c = optimize!(c, x, y, std_noise, mean_θ, std_θ;
			method=method, objective=objective, maxiter=peak_shift_iter,
			regularization=regularization, optimize_mode=optimize_mode, λ=λ,
			verbose=verbose, tol=tol)

		IMs = get_PeakModCP(c, x, mod_peak_num)

		Mod_IMs = optimize!(IMs, x, y, std_noise, peak_mod_mean, peak_mod_std;
					method=bfgs, objective=objective, maxiter=peak_mod_iter,
					regularization=regularization, optimize_mode=optimize_mode,λ=λ,
					verbose=verbose, tol=tol)
		# change_peak_int!.(c.CPs, Mod_IMs[1:end-Int64(have_bg)])
		change_peak_int!.(c.CPs, Mod_IMs)
		# change_c!(c.background, Mod_IMs[end-Int64(have_bg)+1:end])
		c = optimize!(c, x, y, std_noise, mean_θ, std_θ;
			method=method, objective=objective, optimize_mode=optimize_mode,
			maxiter=peak_shift_iter,
			regularization=regularization, λ=λ, verbose=verbose, tol=tol)
	end
	return c
end

function full_optimize!(cp::AbstractVector{<:CrystalPhase}, x::AbstractVector, y::AbstractVector,
						std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
						std_θ::AbstractVector = [1., 1., 5.];
						method::OptimizationMethods, objective::String = "LS",
						optimize_mode::OptimizationMode=Simple,
						regularization::Bool = true,
						loop_num::Int=8,
						peak_shift_iter::Int = 32,
						mod_peak_num::Int = 32,
						peak_mod_mean::AbstractVector = [1.],
						peak_mod_std::AbstractVector = [.5],
						peak_mod_iter::Int=32,
						λ::Float64=1.,
						verbose::Bool = false, tol::Float64 =DEFAULT_TOL)
    pm = PhaseModel(cp)
	pm = full_optimize!(pm, x, y, std_noise, mean_θ, std_θ;
						method=method, objective=objective,
						optimize_mode=Simple,
						regularization=regularization,
						loop_num=loop_num,
						peak_shift_iter=peak_shift_iter,
						mod_peak_num=mod_peak_num,
						peak_mod_mean=peak_mod_mean,
						peak_mod_std=peak_mod_std,
						peak_mod_iter=peak_mod_iter,
						λ=λ,
						verbose=verbose, tol=tol)
	pm.CPs
end

function full_optimize!(cp::CrystalPhase, x::AbstractVector, y::AbstractVector,
	std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
	std_θ::AbstractVector = [1., 1., 5.];
	method::OptimizationMethods, objective::String = "LS",
	optimize_mode::OptimizationMode=Simple,
	regularization::Bool = true,
	loop_num::Int=8,
	peak_shift_iter::Int = 32,
	mod_peak_num::Int = 32,
	peak_mod_mean::AbstractVector = [1.],
	peak_mod_std::AbstractVector = [.5],
	peak_mod_iter::Int=32, λ::Float64=1.,
	verbose::Bool = false, tol::Float64 =DEFAULT_TOL)

	full_optimize!([cp], x, y, std_noise, mean_θ, std_θ;
			method=method, objective=objective,
			optimize_mode=Simple,
			regularization=regularization,
			loop_num=loop_num,
			peak_shift_iter=peak_shift_iter,
			mod_peak_num=mod_peak_num,
			peak_mod_mean=peak_mod_mean,
			peak_mod_std=peak_mod_std,
			peak_mod_iter=peak_mod_iter,
			λ=λ,
			verbose=verbose, tol=tol)
end
"""
    optimize!

    This is the function that each of the node would call on.
    Try to fit a PhaseModel, which comprise a vector of CrystalPhase
	and an optional BackgroundModel to the given spectrum y.

	Return: a PhaseModel object
"""
function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, # Both y and y_uncer will not be further normalized
					std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
					std_θ::AbstractVector = [1., Inf, 5.];
					method::OptimizationMethods=LM, objective::String = "LS",
					optimize_mode::OptimizationMode=Simple,
					maxiter::Int = 32,
					regularization::Bool = true,
					em_loop_num::Integer = 8, λ::Float64=1.,
					verbose::Bool = false, tol::Float64 =DEFAULT_TOL)
	opt_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ,
							maxiter, regularization,
							method, objective, optimize_mode, em_loop_num, λ, verbose, tol)

	optimize!(pm, x, y, y_uncer, opt_stn)
end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector,
				std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
				std_θ::AbstractVector = [1., Inf, 5.];
				method::OptimizationMethods=LM, objective::String = "LS",
				optimize_mode::OptimizationMode=Simple,
				maxiter::Int = 32,
				regularization::Bool = true,
				em_loop_num::Integer = 8, λ::Float64=1.,
				verbose::Bool = false, tol::Float64 =DEFAULT_TOL)
    y_uncer = zero(x)
	optimize!(pm, x, y, y_uncer, std_noise, mean_θ, std_θ,
	          method=method,
			  objective=objective,
	          optimize_mode=optimize_mode,
			  maxiter=maxiter,
	          regularization=regularization,
			  em_loop_num=em_loop_num,
			  verbose=verbose,
			  tol=tol)
end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, opt_stn::OptimizationSettings)
	θ = get_free_params(pm)
	if opt_stn.optimize_mode == Simple
		return simple_optimize!(θ, pm, x, y, y_uncer, opt_stn)
	elseif opt_stn.optimize_mode == EM
		return EM_optimize!(θ, pm, x, y, y_uncer, opt_stn)
    elseif opt_stn.optimize_mode == WithUncer
        return optimize_with_uncertainty!(θ, pm, x, y, opt_stn)
	end
end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings)
	y_uncer = zero(y)
	optimize!(pm, x, y, y_uncer, opt_stn)
end

function simple_optimize!(θ::AbstractVector, pm::PhaseModel,
				   x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, opt_stn::OptimizationSettings)
	if eltype(pm.CPs) <: CrystalPhase && opt_stn.optimize_mode != EM
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
		log_θ = lm_optimize!(log_θ, pm, x, y, y_uncer, opt_stn)
	elseif opt_stn.method == Newton
		log_θ = newton!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == bfgs
		log_θ = BFGS!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == l_bfgs
		log_θ = LBFGS!(log_θ, pm, x, y, opt_stn)
	end

	log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .= @views exp.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
	θ = log_θ
	return reconstruct!(pm, θ)
end

function EM_optimize!(θ::AbstractVector, pm::PhaseModel,
	x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
	opt_stn::OptimizationSettings)

    c = 0 # As existing local to make this thread safe
	std_noise = 0.05
	for i in 1:opt_stn.em_loop_num
		c = simple_optimize!(θ, pm, x, y, y_uncer, opt_stn)
        std_noise = std(y - evaluate(c, get_free_params(c), x))
		θ = get_free_params(c)
		opt_stn = OptimizationSettings{Float64}(opt_stn, std_noise)
		pm = c
	end
	return c
end

function EM_optimize!(θ::AbstractVector, pm::PhaseModel,
	x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings)
    EM_optimize!(θ, pm, x, y, zero(y), opt_stn)
end



function optimize_with_uncertainty!(θ::AbstractVector, pm::PhaseModel,
									x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
									opt_stn::OptimizationSettings)
	# # Background is linear. Hessian is always 0. Need to remove to prevent a weird inexact error

    # t = zeros(Real, length(y) + phase_params)
	# # rr = f(t, phase_log_θ)
	# H = ForwardDiff.hessian(res, phase_log_θ)
	# g = ForwardDiff.gradient(res, phase_log_θ)
	# display(g)
	# display(diag(H))
	# val = res(phase_log_θ)
	# dof = length(x) - length(phase_log_θ)
	# uncer = diag((log((dof-1)*val) - log(dof)) * inverse(H))
	# # uncer = diag((val/dof) * inverse(H))
	# evecs = eigvecs(H)
	# evals = eigvals(H)
    # test_uncer = zeros(Float64, length(evals))
	# for i in 1:size(H, 1)
	# 	for j in 1:size(H, 1)
    #         test_uncer[i] += 1/evals[j] * (evecs[i,j] )^2
	# 	end
	# end
	if eltype(pm.CPs) <: CrystalPhase
		θ = initialize_activation!(θ, pm, x, y)
	end
	# TODO: Don't take log of profile parameters
	θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .= @views log.(θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)]) # tramsform to log space for better conditioning
	log_θ = θ
	(any(isnan, log_θ) || any(isinf, log_θ)) && throw("any(isinf, θ) = $(any(isinf, θ)), any(isnan, θ) = $(any(isnan, θ))")

	# TODO use Match.jl, or just use multiple dispatch on method?
	if opt_stn.method == LM
		log_θ = lm_optimize!(log_θ, pm, x, y, y_uncer, opt_stn)
	elseif opt_stn.method == Newton
		log_θ = newton!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == bfgs
		log_θ = BFGS!(log_θ, pm, x, y, opt_stn)
	elseif opt_stn.method == l_bfgs
		log_θ = LBFGS!(log_θ, pm, x, y, opt_stn)
	end

	# Background is linear. Hessian is always 0. Need to remove to prevent a weird inexact error
	phase_params = get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)
	_, new_bg = reconstruct_BG!(log_θ[phase_params+1:end], pm.background)
	signal = y .- evaluate!(zero(y), new_bg, x)
	phases = PhaseModel(pm.CPs, pm.wildcard, nothing)
	phase_log_θ = log_θ[1:phase_params]

	# This is hessian in log space, TODO: change to real sapce
	if opt_stn.method == LM
		f = get_lm_objective_func(phases, x, signal, y_uncer, opt_stn)
		r = zeros(Real, length(y) + phase_params)
		function res(log_θ)
			sum(abs2, f(r, log_θ))
		end
	else
		res = get_newton_objective_func(pm, x, y, opt_stn)
	end

	# r = zeros(Real, length(y))
	# function t(log_θ)
    #     sum(abs2, y .- evaluate!(r, reconstruct!(pm, exp.(log_θ)), x))
	# end

	H = ForwardDiff.hessian(res, phase_log_θ)
	# val = res(phase_log_θ) * sqrt(2) * opt_stn.priors.std_noise
	val = sum(abs2, y .- evaluate!(zero(x), reconstruct!(pm, exp.(log_θ)), x))
	if opt_stn.verbose
	    println("residual: $(val)")
		display(H)
	end
	# uncer = diag(val / (length(x) - length(phase_log_θ)) * inverse(H))
	uncer = diag(inverse(H))
	log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .= @views exp.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
	θ = log_θ
	pm = reconstruct!(pm, θ)

	# setting fill_angle to zero since structurally determined angles have
	# no uncertainty
	fill_angle = 0
	uncer = get_eight_params(pm.CPs, uncer, fill_angle)
	return pm, uncer
end

function  optimize_with_uncertainty!(θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings)
    optimize_with_uncertainty!(θ, pm, x, y, zero(x), opt_stn)
end


"""
Pass in optimize CrystalPhase arrays and uses Hessian to estimate uncertainty of free parameters
"""
function uncertainty(CPs::AbstractVector{<:CrystalPhase}, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, opt_stn::OptimizationSettings, scaled::Bool=false)
	phase_params = get_param_nums(CPs)
	phase_log_θ = log.(get_free_params(CPs))

	# This is hessian in log space, TODO: change to real sapce
	f = get_lm_objective_func(PhaseModel(CPs, nothing, nothing), x, y, y_uncer, opt_stn)
	r = zeros(Real, length(y) + phase_params)
	function res(log_θ)
		sum(abs2, f(r, log_θ))
	end

	# r = zeros(Real, length(y))
	# function res(log_θ)
	# 	_, new_CP = reconstruct_CPs!(exp.(log_θ), CPs)
    #     sum(abs2, y .- evaluate!(r, new_CP, x))
	# end

	H = ForwardDiff.hessian(res, phase_log_θ)
	l2_res = sum(abs2, y .- evaluate!(zero(x), CPs, x))
	if opt_stn.verbose
	    println("residual: $(val)")
		display(H)
	end


	uncer = scaled ? diag(l2_res / (length(x) - length(phase_log_θ)) * inverse(H)) : diag(inverse(H))
	fill_angle = 0
	uncer = get_eight_params(CPs, uncer, fill_angle)
	return  uncer
end

uncertainty(CPs::AbstractVector{<:CrystalPhase}, x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings, scaled::Bool=false) = uncertainty(CPs, x, y, zero(x), opt_stn, scaled)

function initialize_activation!(θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector)
    new_θ = copy(θ) # make a copy
	start = 1
	for phase in pm.CPs
        param_num = get_param_nums(phase)
		p = evaluate(phase, θ[start:start+param_num-1], x)
		new_θ[start + param_num - 2 - get_param_nums(phase.profile)] = max(0.01, dot(p, y) / sum(abs2, p)) # To avoid crashing with negative value
        start += param_num
	end
	return new_θ
end

function lm_optimize!(log_θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
	                 opt_stn::OptimizationSettings)
	opt_stn.objective == "LS" || error("LM only work with LS for now")

	f = get_lm_objective_func(pm, x, y, y_uncer, opt_stn)
	if opt_stn.regularization
		r = zeros(eltype(log_θ), length(y) + length(log_θ) )
		LM = LevenbergMarquart(f, log_θ, r)
	else
		r = zeros(eltype(log_θ), size(y))
		LM = LevenbergMarquart(f, log_θ, r)
	end

	stn = LevenbergMarquartSettings(min_resnorm = 1e-2, min_res = 1e-3,
						min_decrease = 1e-6, max_iter = opt_stn.maxiter,
						decrease_factor = 7, increase_factor = 10, max_step = .1)

	λ = 1e-6
	OptimizationAlgorithms.optimize!(LM, log_θ, copy(r), stn, λ, Val(opt_stn.verbose))#, false)
	return log_θ
end

function get_lm_objective_func(pm::PhaseModel,
							   x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
							   opt_stn::OptimizationSettings)
	pr = opt_stn.priors
	mean_θ, std_θ = extend_priors(pr, pm)

	function residual!(r::AbstractVector, log_θ::AbstractVector)
		# _sqrt_residual!(pm, log_θ, x, y, r, pr.std_noise)
		_weighted_residual!(pm, log_θ, x, y, y_uncer, r, pr.std_noise)
		# _residual!(pm, log_θ, x, y, r, pr.std_noise)
	end

	function prior!(p::AbstractVector, log_θ::AbstractVector)
		_prior(p, log_θ, mean_θ, std_θ)
	end

	# Regularized cost function
	function f(rp::AbstractVector, log_θ::AbstractVector)
		if (any(isinf, log_θ) || any(isnan, log_θ))
			return Inf
		end
		cp_param_num = get_param_nums(pm.CPs)
		bg_param_num = get_param_nums(pm.background)
		w_param_num = get_param_nums(pm.wildcard)
		θ_cp = log_θ[1:end - w_param_num-bg_param_num]
		θ_w  = log_θ[end - w_param_num-bg_param_num+1 : end-bg_param_num]
		θ_bg = log_θ[end - bg_param_num + 1 : end]
		r = @view rp[1:length(y)] # residual term
		residual!(r, log_θ)
		p = @view rp[length(y)+1:length(y)+cp_param_num] # prior term
		prior!(p, θ_cp)
		wp = @view rp[length(y)+cp_param_num+1:length(y)+cp_param_num+w_param_num]
		lm_prior!(wp, pm.wildcard, θ_w)
		bg_p = @view rp[length(y)+cp_param_num+w_param_num+1:end]
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

	N = LBFGS(get_newton_objective_func(pm, x, y, opt_stn), log_θ, 10, check=false) # default to 10
	N = UnitDirection(N)
	D = DecreasingStep(N, log_θ)
	S = StoppingCriterion(log_θ, dx = tol, rx=tol, maxiter=maxiter, verbose=verbose)
	fixedpoint!(D, log_θ, S)
	return log_θ
end

function BFGS!(log_θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector,
				opt_stn::OptimizationSettings)
	tol, maxiter, verbose = opt_stn.tol, opt_stn.maxiter, opt_stn.verbose

	N = BFGS(get_newton_objective_func(pm, x, y, opt_stn), log_θ, check=false)
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
	λ = opt_stn.λ
	mean_θ, std_θ = extend_priors(pr, pm)
	mean_log_θ = log.(mean_θ)

	function prior(log_θ::AbstractVector)
		bg_param_num = get_param_nums(pm.background)
		w_param_num = get_param_nums(pm.wildcard)
		θ_cp = log_θ[1:end - w_param_num-bg_param_num]
		θ_w  = log_θ[end - w_param_num-bg_param_num+1 : end-bg_param_num]
		θ_bg = log_θ[end - bg_param_num + 1 : end]
		p = zero(eltype(log_θ))
		@inbounds @simd for i in eachindex(θ_cp)
			p += ((θ_cp[i] - mean_log_θ[i]) / (sqrt(2)*std_θ[i]))^2
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
		end_idx = get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)
	    temp_θ = copy(log_θ)
	    @. temp_θ[1:end_idx] = exp(log_θ[1:end_idx])
		if (any(isinf, temp_θ) || any(isnan, temp_θ))
			return Inf
		end
		# @time begin
		r_θ = zeros(promote_type(eltype(log_θ), eltype(x), eltype(y)), length(x))
		evaluate!(r_θ, pm, temp_θ, x)
		# end
		# r_θ = evaluate(pm, temp_θ, x) # reconstruction of phases, TODO: pre-allocate result (one for Dual, one for Float)
		r_θ ./= exp(1) # since we are not normalizing the inputs, this rescaling has the effect that kl(α*y, y) has the optimum at α = 1
		p_θ = prior(log_θ)
		# λ = 1 #TODO: Fix the prior optimization problem and add it to the setting
		# println("p_θ: $(p_θ)")
		# println("kl: $(kl(r_θ, y))")
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
		return (sum(abs2, _prior(p, θ_cp, mean_θ, std_θ))
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

function optimize!(phases::AbstractVector,
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
                   std_θ::AbstractVector = [1., Inf, 5.];
                   method::OptimizationMethods, objective::String = "LS",
				   optimize_mode::OptimizationMode=Simple,
				   maxiter::Int = 32,
				   em_loop_num::Int =1,
				   regularization::Bool = true,
				   λ::Float64 = 1.,
				   verbose::Bool = false, tol::Float64 =DEFAULT_TOL)
	pm = PhaseModel(phases)
	pm = optimize!(pm, x, y, std_noise, mean_θ, std_θ, method=method,
	             objective=objective, maxiter= maxiter,regularization=regularization,
				 optimize_mode=optimize_mode, λ=λ, em_loop_num=em_loop_num,
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


# This actually does not improve much, just cleaner
function _residual!(pm::PhaseModel,
					log_θ::AbstractVector,
					x::AbstractVector, y::AbstractVector,
					r::AbstractVector,
					std_noise::Real)
	end_idx = get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)
	temp_θ = copy(log_θ)
	@. temp_θ[1:end_idx] = exp(log_θ[1:end_idx])
	if (any(isinf, temp_θ) || any(isnan, temp_θ))
		return Inf
	end
	@. r = y
	evaluate_residual!(pm, temp_θ, x, r) # Avoid allocation, put everything in here??
	@. r /= sqrt(2) * std_noise # trade-off between prior and
	# actual residual
	return r
end

function _weighted_residual!(pm::PhaseModel,
					log_θ::AbstractVector,
					x::AbstractVector, y::AbstractVector,
                    y_uncer::AbstractVector,
					r::AbstractVector,
					std_noise::Real)

    end_idx = get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)
    temp_θ = copy(log_θ)
    @. temp_θ[1:end_idx] = exp(log_θ[1:end_idx])
    if (any(isinf, temp_θ) || any(isnan, temp_θ))
		return Inf
    end
    @. r = y
    evaluate_residual!(pm, temp_θ, x, r) # Avoid allocation, put everything in here??
    @. r /= sqrt(2) * sqrt(y_uncer^2 + std_noise^2) # trade-off between prior and
    # actual residual
    return r
end

# using Plots

function _sqrt_residual!(pm::PhaseModel,
						log_θ::AbstractVector,
						x::AbstractVector, y::AbstractVector,
						r::AbstractVector,
						std_noise::Real)

	end_idx = get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)
	temp_θ = copy(log_θ)
	@. temp_θ[1:end_idx] = exp(log_θ[1:end_idx])

	if (any(isinf, temp_θ) || any(isnan, temp_θ))
		return Inf
	end

	@. r = sqrt(y)
	# @. r = y
	evaluate_residual_in_sqrt!(pm, temp_θ, x, r) # Avoid allocation, put everything in here??
	# if eltype(r) <: Float64
	# 	plt = plot(x, sqrt.(y))
	# 	plot!(x, r)
	# 	display(plt)
	# end
	r ./= sqrt(2) * std_noise # trade-off between prior and

	# actual residual
	return r
end

########################### parameter helpers ##################################
function check_objective(objective::String)
	objective in ALLOWED_OBJECTIVE || error("objective $(objective) not a allowed objective string")
end