"""
    `optimize!`

    Function for optimizing `PhaseModel`s.
	There are multiple functions that handles more primitive object, e.g. Vector{CrystalPhase}.
	Returns a optimized `PhaseModel`.
"""
function optimize! end

"""
   `full_optimize!`

    Function for optimizing `PhaseModel`s. This function allows change in peak height ratios.
	There are multiple functions that handles more primitive object, e.g. Vector{CrystalPhase}.
	Returns a optimized `PhaseModel`.
"""
function full_optimize! end


function fit_phases(phases::AbstractVector{<:CrystalPhase}, x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean_θ::AbstractVector = [1. , .5, .2], std_θ::AbstractVector = [.005, .5, .2];
                   maxiter::Int = 5122, regularization::Bool = true)
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

function fit_amorphous(W::Wildcard, BG::Background, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
					std_noise::Real;
					method::OptimizationMethod,
					objective::String = "LS",
					optimize_mode::OptimizationMode=Simple,
					maxiter::Int = 512,
					regularization::Bool = true,
					em_loop_num::Integer = 8,
					λ::Float64 = 1.,
					verbose::Bool = false, tol::Float64 = DEFAULT_TOL)

    pm = PhaseModel(W, BG)
	opt_stn = OptimizationSettings{Float64}(std_noise, [1., .5, .2], [.005, .5, .2],
											maxiter, regularization,
											method, objective, optimize_mode, em_loop_num, λ, verbose, tol)

	y ./= maximum(y) * 2
	opt_pm = optimize!(pm, x, y, y_uncer, opt_stn, opt_stn.optimize_mode)
	return opt_pm
end

# TODO: allow change in some parameters
"""
    full_optimize!(pm::PhaseModel, x, y, std_noise, mean_θ, std_θ; kwargs...)

Optimize a `PhaseModel` while also allowing the peak heights to change. Each of the
`loop_num` loops refines the lattice, activation, width and background (`optimize!`),
then fits multiplicative height factors for the first `mod_peak_num` peaks of each
phase with the positions fixed, and refines the lattice again. The phases passed in
are not modified; the fitted peak heights are in the returned model.

The height factors are relative to the reference intensities of the phases passed in
(`pm.CPs`): every loop re-solves them from the reference intensities, so the log-normal
prior (`peak_mod_mean`, `peak_mod_std`) bounds the total deviation from the reference,
however many loops are run.
"""
function full_optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector,
						std_noise::Real, mean_θ::AbstractVector = [1., .5, .2], std_θ::AbstractVector = [.005, .5, .2];
						method::OptimizationMethod = LM,
						optimize_mode::OptimizationMode=Simple,
						objective::String = "LS",
						regularization::Bool = true,
						loop_num::Int=8,
						peak_shift_iter::Int = 32,
						mod_peak_num::Int = 32,
						peak_mod_mean::AbstractVector = [1.],
						peak_mod_std::AbstractVector = [.5],
						peak_mod_iter::Int=32,
						analytic_peak_mod::Bool = true, # solve peak heights with LinearPeakMod (no autodiff)
						λ::Float64=1.,
						verbose::Bool = false, tol::Float64 =DEFAULT_TOL)

	# have_bg = !isnothing(pm.background)
	# Work on copies of the peaks so change_peak_int! does not modify the caller's phases
	ref_peaks = [CP.peaks for CP in pm.CPs]
	c = PhaseModel(copy_peaks.(pm.CPs), pm.wildcard, pm.background)
	for i in 1:loop_num
		c = optimize!(c, x, y, std_noise, mean_θ, std_θ;
			method=method, objective=objective, maxiter=peak_shift_iter,
			regularization=regularization, optimize_mode=optimize_mode, λ=λ,
			verbose=verbose, tol=tol)

		# Fit the height factors relative to the reference intensities, so the prior
		# applies to the total deviation instead of compounding over loops
		c = PhaseModel(with_peaks.(c.CPs, ref_peaks), c.wildcard, c.background)
		IMs = get_PeakModCP(c, x, mod_peak_num)

		if analytic_peak_mod && objective == "LS" && optimize_mode == Simple
			Mod_IMs = optimize!(LinearPeakMod(IMs, y, peak_mod_mean, peak_mod_std), std_noise;
						maxiter=peak_mod_iter, verbose=verbose, tol=tol)
		else
			Mod_IMs = optimize!(IMs, x, y, std_noise, peak_mod_mean, peak_mod_std;
						method=bfgs, objective=objective, maxiter=peak_mod_iter,
						regularization=regularization, optimize_mode=optimize_mode,λ=λ,
						verbose=verbose, tol=tol)
		end
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
						method::OptimizationMethod, objective::String = "LS",
						optimize_mode::OptimizationMode=Simple,
						regularization::Bool = true,
						loop_num::Int=8,
						peak_shift_iter::Int = 32,
						mod_peak_num::Int = 32,
						peak_mod_mean::AbstractVector = [1.],
						peak_mod_std::AbstractVector = [.5],
						peak_mod_iter::Int=32,
						analytic_peak_mod::Bool = true,
						λ::Float64=1.,
						verbose::Bool = false, tol::Float64 =DEFAULT_TOL)
    pm = PhaseModel(cp)
	pm = full_optimize!(pm, x, y, std_noise, mean_θ, std_θ;
						method=method, objective=objective,
						optimize_mode=optimize_mode,
						regularization=regularization,
						loop_num=loop_num,
						peak_shift_iter=peak_shift_iter,
						mod_peak_num=mod_peak_num,
						peak_mod_mean=peak_mod_mean,
						peak_mod_std=peak_mod_std,
						peak_mod_iter=peak_mod_iter,
						analytic_peak_mod=analytic_peak_mod,
						λ=λ,
						verbose=verbose, tol=tol)
	pm.CPs
end

function full_optimize!(cp::CrystalPhase, x::AbstractVector, y::AbstractVector,
	std_noise::Real, mean_θ::AbstractVector = [1., 1., .2],
	std_θ::AbstractVector = [1., 1., 5.];
	method::OptimizationMethod, objective::String = "LS",
	optimize_mode::OptimizationMode=Simple,
	regularization::Bool = true,
	loop_num::Int=8,
	peak_shift_iter::Int = 32,
	mod_peak_num::Int = 32,
	peak_mod_mean::AbstractVector = [1.],
	peak_mod_std::AbstractVector = [.5],
	peak_mod_iter::Int=32, λ::Float64=1.,
	analytic_peak_mod::Bool = true,
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
			analytic_peak_mod=analytic_peak_mod,
			λ=λ,
			verbose=verbose, tol=tol)
end

"""
    optimize!(P::LinearPeakMod, std_noise; maxiter, tol, verbose)

Optimize the peak-height factors w = exp(u) of a `LinearPeakMod` with damped
Gauss-Newton (Levenberg-Marquardt) using the analytic gradient and Gauss-Newton
Hessian. Minimizes the same objective as the BFGS/`PeakModCP` route with "LS":

    F(u) = ‖y - c - B w‖² / (2 std_noise²) + Σ ((u - mean_log_θ) / (√2 std_θ))²

Returns a vector of `PeakModCP`s holding the optimized height factors, one per phase.

Notes:
- Like the BFGS route, the prior is always applied (there is no `regularization`
  switch) and measurement uncertainty `y_uncer` is not used.
- `tol` is a stopping threshold on the objective decrease of an accepted step,
  relative to max(1, F), or on the largest step component in u; this differs from
  the `dx`/`rx` criteria of the BFGS route.
- `full_optimize!` uses this route when `analytic_peak_mod = true` (default),
  `objective == "LS"` and `optimize_mode == Simple`, and falls back to BFGS otherwise.
"""
function optimize!(P::LinearPeakMod, std_noise::Real;
				   maxiter::Int = 32, tol::Real = DEFAULT_TOL, verbose::Bool = false)
	σ² = std_noise^2
	s² = P.std_θ .^ 2
	u = log.(reduce(vcat, [get_free_params(IM) for IM in P.IMs]))

	function objective(u)
		w = exp.(u)
		(P.rtr - 2dot(w, P.Btr) + dot(w, P.BtB, w)) / (2σ²) + sum((u .- P.mean_log_θ).^2 ./ (2s²))
	end

	F = objective(u)
	damping = 1e-6
	for i in 1:maxiter
		w = exp.(u)
		g = w .* (P.BtB * w .- P.Btr) ./ σ² .+ (u .- P.mean_log_θ) ./ s²
		H = (w * w') .* P.BtB ./ σ² + Diagonal(1 ./ s²) # Gauss-Newton Hessian, positive definite
		accepted = false
		while damping < 1e10
			δ = -(Symmetric(H + damping * Diagonal(diag(H))) \ g)
			F_new = objective(u .+ δ)
			if F_new < F
				u .+= δ
				decrease = F - F_new
				F = F_new
				damping = max(damping / 7, 1e-12)
				accepted = true
				verbose && println("LinearPeakMod iter $i: objective = $F")
				(decrease <= tol * max(1, F) || maximum(abs, δ) <= tol) && return reconstruct_IMs(P, exp.(u))
				break
			end
			damping *= 10
		end
		accepted || break # no decrease possible, at a minimum up to numerical precision
	end
	reconstruct_IMs(P, exp.(u))
end




function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, # Both y and y_uncer will not be further normalized
					std_noise::Real, mean_θ::AbstractVector = [1., .5, .2], std_θ::AbstractVector = [.005, .5, .2];
					method::OptimizationMethod=LM,
					objective::String = "LS",
					optimize_mode::OptimizationMode=Simple,
					maxiter::Int = 512,
					regularization::Bool = true,
					em_loop_num::Integer = 8, λ::Float64=1.,
					verbose::Bool = false, tol::Float64 =DEFAULT_TOL)
	opt_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ,
							maxiter, regularization,
							method, objective, optimize_mode, em_loop_num, λ, verbose, tol)

	optimize!(pm, x, y, y_uncer, opt_stn, opt_stn.optimize_mode)
end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector,
				std_noise::Real, mean_θ::AbstractVector = [1., .5, .2], std_θ::AbstractVector = [.005, .5, .2];
				method::OptimizationMethod=LM, objective::String = "LS",
				optimize_mode::OptimizationMode=Simple,
				maxiter::Int = 512,
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

# function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, opt_stn::OptimizationSettings)
# 	θ = get_free_params(pm)
# 	if opt_stn.optimize_mode == Simple
# 		return simple_optimize!(θ, pm, x, y, y_uncer, opt_stn)
# 	elseif opt_stn.optimize_mode == EM
# 		return EM_optimize!(θ, pm, x, y, y_uncer, opt_stn)
#     elseif opt_stn.optimize_mode == WithUncer
#         return optimize_with_uncertainty!(θ, pm, x, y, opt_stn)
# 	end
# end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, opt_stn::OptimizationSettings)
	optimize!(pm, x, y, y_uncer, opt_stn, opt_stn.optimize_mode) # Dispatch on optimize_mode, which gives different output that affects type stability
end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, opt_stn::OptimizationSettings, mode::_Simple)
    θ = get_free_params(pm)
    simple_optimize!(θ, pm, x, y, y_uncer, opt_stn)
end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, opt_stn::OptimizationSettings, mode::_EM)
    θ = get_free_params(pm)
    EM_optimize!(θ, pm, x, y, y_uncer, opt_stn)
end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, opt_stn::OptimizationSettings, mode::_WithUncer)
    θ = get_free_params(pm)
    optimize_with_uncertainty!(θ, pm, x, y, opt_stn)
end

function optimize!(pm::PhaseModel, x::AbstractVector, y::AbstractVector, opt_stn::OptimizationSettings)
	y_uncer = zero(y)
	optimize!(pm, x, y, y_uncer, opt_stn)
end

function simple_optimize!(θ::AbstractVector, pm::PhaseModel,
				   x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, opt_stn::OptimizationSettings)
	# TODO: Add option for whether to estimate this, otherwise it will casue issue with loop methods (EM, full_optimize)
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
	elseif opt_stn.method == dogleg
		log_θ = dogleg_optimize!(log_θ, pm, x, y, y_uncer, opt_stn)
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

# FIXME: Really high allocation counts
function EM_optimize!(θ::AbstractVector, pm::PhaseModel,
	x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
	opt_stn::OptimizationSettings)

    c = 0 # As existing local to make this thread safe
	std_noise = 0.05
	xrd_temp = zero(y)

	for i in 1:opt_stn.em_loop_num
		c = simple_optimize!(θ, pm, x, y, y_uncer, opt_stn)

		if i != opt_stn.em_loop_num
		    evaluate!(xrd_temp, c, get_free_params(c), x)
		    std_noise = std(y .- xrd_temp)
            xrd_temp .= 0
		    opt_stn = OptimizationSettings{Float64}(opt_stn, std_noise, round(Int64, 32)) # Decreasing iters
		    pm = c
		end
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
	elseif opt_stn.method == dogleg
		log_θ = dogleg_optimize!(log_θ, pm, x, y, y_uncer, opt_stn)
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
	if opt_stn.method in (LM, dogleg) # least-squares objective
		f = get_lm_objective_func(phases, x, signal, y_uncer, opt_stn)
		r = zeros(Real, length(y) + phase_params)
		function res(log_θ)
			sum(abs2, f(r, log_θ))
		end
	else
		res = get_newton_objective_func(pm, x, y, opt_stn)
	end

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
`uncertainty`

Pass in optimize CrystalPhase arrays and uses Hessian to estimate uncertainty of free parameters.
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
	temp = zero(x)
	for phase in pm.CPs
        param_num = get_param_nums(phase)
		evaluate!(temp, phase, θ[start:start+param_num-1], x)
		# Add maximum of 1.0, because it should not be too much larger than that. Most likley becuase of strong background
		new_θ[start + param_num - 2 - get_param_nums(phase.profile)] = min(1.0, max(0.01, dot(temp, y) / sum(abs2, temp))) # To avoid crashing with negative value
		temp .= 0
        start += param_num
	end
	return new_θ
end

function lm_optimize!(log_θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
	                 opt_stn::OptimizationSettings)
	opt_stn.objective == "LS" || error("LM only work with LS for now")

	f = get_lm_objective_func(pm, x, y, y_uncer, opt_stn)
	r = zeros(eltype(log_θ), opt_stn.regularization ? length(y) + length(log_θ) : length(y))

	stn = LevenbergMarquartSettings(min_resnorm = 1e-2, min_res = 1e-3,
						min_decrease = 1e-6, max_iter = opt_stn.maxiter,
						decrease_factor = 7, increase_factor = 10, max_step = .1)

	λ = 1e-6
	np = length(log_θ) - get_param_nums(pm.background)
	if is_linear(pm.background) && np > 0
		# background is linear: AD only over phase/wildcard params, constant Jacobian for the rest
		J_bg = background_jacobian(pm.background, x, y, y_uncer, length(r), opt_stn)
		LM = SplitJacobianLM(f, log_θ, r, np, J_bg)
		lm_loop!(LM, log_θ, copy(r), stn, λ, Val(opt_stn.verbose))
	else
		LM = LevenbergMarquart(f, log_θ, r)
		OptimizationAlgorithms.optimize!(LM, log_θ, copy(r), stn, λ, Val(opt_stn.verbose))#, false)
	end
	return log_θ
end

# Largest relative change of a lattice parameter (lengths and angles) allowed in one
# dogleg_optimize! call, as a box bound in log space
const DOGLEG_MAX_STRAIN = 0.05

# Same least-squares problem as lm_optimize! (residual, priors, log-space parameters),
# solved with the Dogleg trust-region method of LeastSquaresOptim. In the benchmarks in
# benchmark/FINDINGS.md it took far fewer iterations than the LM above and recovered
# lattice parameters at least as well. The Jacobian comes from ForwardDiff over all
# parameters, including linear backgrounds.
# Unlike the LM, whose steps are capped at 0.1 in log space, the trust region can take
# large steps, which let wrong phases strain far to fit part of a pattern during tree
# search; lattice parameters are therefore bounded to ±DOGLEG_MAX_STRAIN of the start.
function dogleg_optimize!(log_θ::AbstractVector, pm::PhaseModel, x::AbstractVector, y::AbstractVector,
                          y_uncer::AbstractVector, opt_stn::OptimizationSettings)
	opt_stn.objective == "LS" || error("dogleg only works with LS for now")

	f = get_lm_objective_func(pm, x, y, y_uncer, opt_stn)
	n_res = opt_stn.regularization ? length(y) + length(log_θ) : length(y)
	function f!(r, θ)
		out = f(r, θ)
		out isa Number && fill!(r, Inf) # the objective returns a scalar Inf for non-finite parameters
		r
	end
	lower, upper = fill(-Inf, length(log_θ)), fill(Inf, length(log_θ))
	if !isnothing(pm.CPs) && eltype(pm.CPs) <: CrystalPhase
		start = 1
		for cp in pm.CPs
			lattice = start:start+cp.cl.free_param-1
			lower[lattice] .= log_θ[lattice] .+ log(1 - DOGLEG_MAX_STRAIN)
			upper[lattice] .= log_θ[lattice] .+ log(1 + DOGLEG_MAX_STRAIN)
			start += get_param_nums(cp)
		end
	end
	θ0 = copy(log_θ)
	problem = LeastSquaresOptim.LeastSquaresProblem(x = log_θ, f! = f!, output_length = n_res,
	                                                autodiff = :forward)
	result = LeastSquaresOptim.optimize!(problem, LeastSquaresOptim.Dogleg(); lower, upper,
	                                     iterations = opt_stn.maxiter, show_trace = opt_stn.verbose)
	log_θ .= all(isfinite, result.minimizer) ? result.minimizer : θ0
	return log_θ
end

# Constant Jacobian of the LM residual vector (see get_lm_objective_func) w.r.t.
# the coefficients of a linear background: data rows and background prior rows
function background_jacobian(B::AbstractBackground, x::AbstractVector, y::AbstractVector,
                             y_uncer::AbstractVector, n_rows::Int, opt_stn::OptimizationSettings)
	nbg = get_param_nums(B)
	w = @. 1 / (sqrt(2) * sqrt(y_uncer^2 + opt_stn.priors.std_noise^2)) # as in _weighted_residual!
	J = zeros(n_rows, nbg)
	J[1:length(y), :] .= -w .* background_basis(B, x)
	if opt_stn.regularization # background prior p = Λ c
		Λ = zeros(nbg)
		lm_prior!(Λ, B, ones(nbg))
		J[end-nbg+1:end, :] .= Diagonal(Λ)
	end
	J
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

	N = BFGS(get_newton_objective_func(pm, x, y, opt_stn), log_θ, check=true)
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
                   method::OptimizationMethod, objective::String = "LS",
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
					method::OptimizationMethod, objective::String = "LS",
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
