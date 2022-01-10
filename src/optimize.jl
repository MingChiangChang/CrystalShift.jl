function fit_phases(phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real = .1, mean_θ::AbstractVector = [1. ,.2],
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
	# Preallocate or parrellel
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
                   std_noise::Real = .001, mean_θ::AbstractVector = [1., .2],
                   std_θ::AbstractVector = [1., 5.];
                   maxiter::Int = 32, regularization::Bool = true)
    θ = get_parameters(phases)

	if length(mean_θ) == 3 #2 
		# Different prior for different crystals?
	    mean_θ, std_θ = extend_priors(mean_θ, std_θ, phases)
	end
	println("mean_θ = $(mean_θ)")
	println("std_θ = $(std_θ)")

	length_check(phases, mean_θ, std_θ) || error("number of parameter must match number of terms in the prior")

    θ = optimize!(θ, phases, x, y, std_noise, mean_θ, std_θ,
	          maxiter = maxiter, regularization = regularization)

    for (i, cp) in enumerate(phases)
        phases[i] = CrystalPhase(cp, θ)
		deleteat!(θ, collect(1:get_param_nums(phases[i])))
	end
	return phases
end

function length_check(phases::AbstractVector, mean_θ::AbstractVector, std_θ::AbstractVector)
	# println("sum: $(sum([phase.cl.free_param + 1 for phase in phases]))")
	# println("")
    (sum([phase.cl.free_param + 2 for phase in phases]) == length(mean_θ) == length(std_θ) )
end

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

function extend_priors(mean_θ::AbstractVector, std_θ::AbstractVector,
	                    phases::AbstractVector{<:CrystalPhase})
	extend_priors(mean_θ, std_θ, [phase.cl for phase in phases])
end

function extend_priors(mean_θ::AbstractVector, std_θ::AbstractVector,
	                    phases::AbstractVector{<:Crystal})
	totl_params = sum([(phase.free_param + 2) for phase in phases])
	full_mean_θ = zeros(totl_params)
	full_std_θ = zeros(totl_params)
	start = 1
	for phase in phases
		n = phase.free_param
		full_mean_θ[start:start+n-1] = mean_θ[1].*get_free_params(phase)
		full_std_θ[start:start+n-1] = std_θ[1].*get_free_params(phase)#repeat(std_θ[1, :], n)
		full_mean_θ[start + n:start + n + 1] = mean_θ[2:3]
		full_std_θ[start + n:start + n + 1] = std_θ[2:3]
		start += (n+2)
    end
	return full_mean_θ, full_std_θ
end

# Single phase situation. Put phase into [phase].
function optimize!(phase::CrystalPhase, x::AbstractVector, y::AbstractVector,
                   std_noise::Real = .01, mean_θ::AbstractVector = [1., .2],
                   std_θ::AbstractVector = [.1, 1.];
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

# TODO create prior for each crystal phases
# TODO rename variables that are in log space

function optimize!(θ::AbstractVector, phases::AbstractVector{<:CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean_θ::AbstractVector,
                   std_θ::AbstractVector; maxiter::Int = 32,
                   regularization::Bool = true)
	# params = copy(θ)
    function residual!(r::AbstractVector, log_θ::AbstractVector)
        params = exp.(log_θ) # make a copy

        @. r = y
		res!(phases, params, x, r) # Avoid allocation, put everything in here??
		# r -= reconstruct!((phases,), (params,), x, temp)
		r ./= sqrt(2) * std_noise # trade-off between prior and
		                          # actual residual
        return r
	end

	# The lower symmetry phases have better fitting power and thus
	# should be punished more by the prior
    function prior!(p::AbstractVector, log_θ::AbstractVector)
		# θ_c = remove_act_from_θ(log_θ, phases)
		μ = log.(mean_θ)
		# @. p = (θ_c - μ) / (sqrt(2)*std_θ)
		@. p = (log_θ - μ) / (sqrt(2)*std_θ)
	end

	# Regularized cost function
    function f(rp::AbstractVector, log_θ::AbstractVector)
    	r = @view rp[1:length(y)] # residual term
    	residual!(r, log_θ)
    	p = @view rp[length(y)+1:end] # prior term
    	prior!(p, log_θ)
		#println("Norm: $(norm(r)) Prior:$(norm(p))")
    	return rp
    end

    θ = initialize_activation!(θ, phases, x, y)

    @. θ = log(θ) # tramsform to log space for better conditioning
	log_θ = θ
    (any(isnan, log_θ) || any(isinf, log_θ)) && throw("any(isinf, θ) = $(any(isinf, θ)), any(isnan, θ) = $(any(isnan, θ))")
    println(log_θ)
    if regularization
		r = zeros(eltype(log_θ), length(y) + length(log_θ) ) # Reason?? - size(phases, 1)
        LM = LevenbergMarquart(f, log_θ, r)
	else
		r = zeros(eltype(log_θ), size(y))
        LM = LevenbergMarquart(residual!, log_θ, r)
	end
	
    stn = LevenbergMarquartSettings(min_resnorm = 1e-2, min_res = 1e-3,
						min_decrease = 1e-8, max_iter = maxiter,
						decrease_factor = 7, increase_factor = 10, max_step = 0.1)
	λ = 1e-6

	OptimizationAlgorithms.optimize!(LM, log_θ, copy(r), stn, λ, Val(false))
	@. θ = exp(log_θ) # transform back
	return θ
end
