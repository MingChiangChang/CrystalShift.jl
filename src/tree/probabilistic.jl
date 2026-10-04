# computes log marginal likelihood of θ given (x, y) based on the Laplace approximation
# NOTE: the input θ should be an (approximate) local minimum with respect to θ
# mean_θ, std_θ are the mean and standard deviation of the prior Gaussian distribution of θ
# std_noise is the standard deviation of the noise, can be real number or vector
# if std_noise is a vector, each element corresponds the std of an element in y
# TODO: Clean up repeated code and use dispatch instead
function get_probabilities(results::AbstractVector{<:Node},
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        y_uncer::AbstractVector,
        std_noise::Real,
        mean_θ::AbstractVector{<:Real},
        std_θ::AbstractVector{<:Real};
        objective::AbstractObjective = LeastSquares(),
        renormalize::Bool = true,
        normalization_constant::Real = 1.,
        n_sigma::Real = 1.)

    neg_log_prob = zeros(length(results))
    std_noise = n_sigma * std_noise

    @threads for i in eachindex(results) # Phiddle not happy with using @threads here
        θ = get_free_params(results[i].phase_model)
        full_mean_θ, full_std_θ = extend_priors(mean_θ, std_θ, results[i].phase_model.CPs)
        neg_log_prob[i] = approximate_negative_log_evidence(results[i], θ, x, y, y_uncer,
                                                            std_noise, full_mean_θ, full_std_θ, objective)
    end

    if renormalize
        neg_log_prob ./= minimum(neg_log_prob) * std_noise # Renormalize
        neg_log_prob .*= normalization_constant
    end

    log_normalization = logsumexp(-neg_log_prob)  # numerically stable computation

    return @. exp(-(neg_log_prob + log_normalization)) # i.e. prob / sum(prob)
end

function get_probabilities(results::AbstractVector{<:Node},
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        std_noise::Real,
        mean_θ::AbstractVector{<:Real},
        std_θ::AbstractVector{<:Real};
        objective::AbstractObjective = LeastSquares(),
        renormalize::Bool = true,
        normalization_constant::Real = 1.,
        n_sigma::Real = 1.)

    y_uncer = zero(y)
    get_probabilities(results, x, y, y_uncer, std_noise, mean_θ, std_θ,
                      objective=objective,
                      renormalize=renormalize,
                      normalization_constant=normalization_constant,
                      n_sigma=n_sigma)
end



function get_probabilities(results::AbstractVector{<:Node},
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        y_uncer::AbstractVector,
        mean_θ::AbstractVector{<:Real},
        std_θ::AbstractVector{<:Real};
        objective::AbstractObjective = LeastSquares(),
        renormalize::Bool = true,
        normalization_constant::Real = 1.,
        n_sigma::Real = 1.)

    std_noise = n_sigma * minimum([std(y.- evaluate!(zero(x), results[i].phase_model, x)) for i in eachindex(results)])
    get_probabilities(results, x, y, y_uncer, std_noise, mean_θ, std_θ,
                      objective=objective,
                      renormalize=renormalize,
                      normalization_constant=normalization_constant)
end

function get_probabilities(results::AbstractVector{<:Node},
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        mean_θ::AbstractVector{<:Real},
        std_θ::AbstractVector{<:Real};
        objective::AbstractObjective = LeastSquares(),
        renormalize::Bool = true,
        normalization_constant::Real = 1.,
        n_sigma::Real = 1.)
    y_uncer = zero(y)
    std_noise = n_sigma * minimum([std(y.- evaluate!(zero(x), results[i].phase_model, x)) for i in eachindex(results)])
    get_probabilities(results, x, y, y_uncer, std_noise, mean_θ, std_θ,
                      objective=objective,
                      renormalize=renormalize,
                      normalization_constant=normalization_constant)
end


# λ controls the weight of the least-squares regularization term for kl objective
function approximate_negative_log_evidence(node::Node, θ::AbstractVector, x::AbstractVector,
        y::AbstractVector, std_noise::Real,
        mean_θ::RealOrVec, std_θ::RealOrVec,
        objective::AbstractObjective, λ::Real = 1e-6,
        verbose::Bool = false)
    mean_log_θ = log.(mean_θ)
    y_uncer = zero(y)

    f = if objective isa LeastSquares
        function (log_θ)
            ls_objective(node.phase_model, log_θ, x, y, y_uncer, std_noise, mean_log_θ, std_θ)
        end
    elseif objective isa KullbackLeilber
        function (log_θ)
            kl_objective(node.phase_model, log_θ, x, y, mean_log_θ, std_θ, λ)
        end
    end

    # tramsform to log space for better conditioning
    θ[1:get_param_nums(node.phase_model.CPs)]= log.(θ[1:get_param_nums(node.phase_model.CPs)])
    log_θ = θ
    # IDEA: chek norm of gradient. If it exceeds a threshold, apply fine-tuning.
    # newton!(f, log_θ)
    return approximate_negative_log_evidence(f, log_θ, verbose)
end

function approximate_negative_log_evidence(node::Node,
        θ::AbstractVector,
        x::AbstractVector,
        y::AbstractVector,
        y_uncer::AbstractVector,
        std_noise::Real,
        mean_θ::RealOrVec,
        std_θ::RealOrVec,
        objective::AbstractObjective,
        λ::Real = 1e-6,
        verbose::Bool = false)

    mean_log_θ = log.(mean_θ)

    f = if objective isa LeastSquares
        function (log_θ)
            ls_objective(node.phase_model, log_θ, x, y, y_uncer, std_noise, mean_log_θ, std_θ)
        end
    elseif objective isa KullbackLeilber
        function (log_θ)
            kl_objective(node.phase_model, log_θ, x, y, mean_log_θ, std_θ, λ)
        end
    end

    # tramsform to log space for better conditioning
    θ[1:get_param_nums(node.phase_model.CPs)]= log.(θ[1:get_param_nums(node.phase_model.CPs)])
    log_θ = θ
    # IDEA: chek norm of gradient. If it exceeds a threshold, apply fine-tuning.
    # newton!(f, log_θ)
    return approximate_negative_log_evidence(f, log_θ, verbose)
end

# function get_probabilities(results::AbstractVector{<:Node},
#         x::AbstractVector{<:Real},
#         y::AbstractVector{<:Real},
#         mean_θ::AbstractVector{<:Real},
#         std_θ::AbstractVector{<:Real};
#         objective::AbstractObjective = LeastSquares(),
#         renormalize::Bool = true,
#         normalization_constant::Real = 1.0)
# 
#     neg_log_prob = zeros(length(results))
#     std_noise = minimum([std(y.- evaluate!(zero(x), results[i].phase_model, x)) for i in eachindex(results)])
#     # println(std_noise)
# 
#     for i in 1:length(results)
#         # θ = get_free_params_w_std_noise(results[i].phase_model, x, y)
#         θ = get_free_params(results[i].phase_model)
#         #θ[end] = std_noise
#         full_mean_θ, full_std_θ = extend_priors(mean_θ, std_θ, results[i].phase_model.CPs)
#         neg_log_prob[i] = approximate_negative_log_evidence(results[i], θ, x, y,
#                                                             std_noise, full_mean_θ, full_std_θ, objective)
#     end
# 
#     # for i in eachindex(results)
#     # 	plt = plot(x, y)
#     # 	plot!(x, evaluate!(zero(x), results[i].phase_model, x), title="$(neg_log_prob[i])")
#     # 	display(plt)
#     # end
#     # println(minimum(neg_log_prob))
#     if renormalize
#         # renormalize
#         neg_log_prob ./= minimum(neg_log_prob) * std_noise# * 0.07 #  This makes good calibration curve but pushes things to both sides
#         # neg_log_prob .*= std_noise
#         # neg_log_prob .-= minimum(neg_log_prob)
#         neg_log_prob .*= normalization_constant
#     end
#     log_normalization = logsumexp(-neg_log_prob)  # numerically stable computation
#     return @. exp(-(neg_log_prob + log_normalization)) # i.e. prob / sum(prob)
# end

function approximate_negative_log_evidence(node::Node, θ::AbstractVector, x::AbstractVector,
        y::AbstractVector,
        mean_θ::RealOrVec, std_θ::RealOrVec,
        objective::AbstractObjective, λ::Real = 1e-6,
        verbose::Bool = false)
    mean_log_θ = log.(mean_θ)

    f = if objective isa LeastSquares
        function (log_θ)
            ls_objective(node.phase_model, log_θ, x, y, mean_log_θ, std_θ)
        end
    else
        error("unimplemented")
    end

    # tramsform to log space for better conditioning
    θ[1:get_param_nums(node.phase_model.CPs)]= log.(θ[1:get_param_nums(node.phase_model.CPs)])
    log_θ = θ

    return approximate_negative_log_evidence(f, log_θ, verbose)
end

function get_free_params_w_std_noise(phases::AbstractVector{<:AbstractPhase}, x::AbstractVector, y::AbstractVector)
    θ = zeros(promote_type(eltype(x), eltype(y)), get_param_nums(phases) + 1)
    θ[1:end-1] .= get_free_params(phases)
    θ[end] = std(y - evaluate!(zero(x), phases, x)) # std_noise
    θ
end

get_free_params_w_std_noise(node::Node, x::AbstractVector, y::AbstractVector) = get_free_params_w_std_noise(node.phase_model.CPs, x, y)
get_free_params_w_std_noise(pm::PhaseModel, x::AbstractVector, y::AbstractVector) = get_free_params_w_std_noise(pm.CPs, x, y)

function ls_objective(pm::PhaseModel, log_θ::AbstractVector, x::AbstractVector, y::AbstractVector,
        mean_log_θ::AbstractVector, std_θ::AbstractVector)
    ls_residual(pm, log_θ, x, y) + ls_regularizer(pm, log_θ[1:end-1], mean_log_θ, std_θ)
end

function ls_residual(pm::PhaseModel, log_θ::AbstractVector, x::AbstractVector, y::AbstractVector)
    r = zeros(promote_type(eltype(log_θ), eltype(x), eltype(y)), length(x))
    # r = _residual!(PM, log_θ, x, y, r, std_noise)
    log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .= @views exp.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
    if (any(isinf, log_θ) || any(isnan, log_θ))
        log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .=  @views log.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
        return Inf
    end
    @. r = y
    evaluate_residual!(pm, log_θ, x, r) # Avoid allocation, put everything in here??
    r ./= sqrt(2) * log_θ[end] # trade-off between prior and

    # actual residual
    log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .=  @views log.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
    return sum(abs2, r)
end

# function ls_objective(pm::PhaseModel, log_θ::AbstractVector, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
# 					mean_log_θ::AbstractVector, std_θ::AbstractVector)
# 	ls_residual(pm, log_θ, x, y, y_uncer) + ls_regularizer(pm, log_θ[1:end-1], mean_log_θ, std_θ)
# end

# function ls_residual(pm::PhaseModel, log_θ::AbstractVector, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector)
# 	r = zeros(promote_type(eltype(log_θ), eltype(x), eltype(y)), length(x))
# 	# r = _residual!(PM, log_θ, x, y, r, std_noise)
# 	log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .= @views exp.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
# 	if (any(isinf, log_θ) || any(isnan, log_θ))
# 		log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .=  @views log.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
# 		return Inf
# 	end
# 	@. r = y
# 	evaluate_residual!(pm, log_θ, x, r) # Avoid allocation, put everything in here??
# 	@. r /= sqrt(2) * sqrt(y_uncer^2 + log_θ[end]^2) # trade-off between prior and

# 	# actual residual
# 	log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)] .=  @views log.(log_θ[1:get_param_nums(pm.CPs)+get_param_nums(pm.wildcard)])
# 	return sum(abs2, r)
# end



# NOTE assumes θ is stationary point of θ (i.e. ∇f = 0)
# computes approximation to \int exp(-f(θ)) dθ via Laplace approximation
function approximate_negative_log_evidence(f, θ, verbose::Bool = false)
    # calculate marginal likelihood
    d = length(θ)
    val = f(θ)
    H = ForwardDiff.hessian(f, θ)
    Σ = inverse(H) # reinterpret Hessian of minimization problem as inverse of covariance matrix
    if verbose
        println("in approximate_negative_log_evidence")
        display(eigvals(Matrix(Σ)))
        display(Matrix(Σ))
    end
    # println("val: $(val)")
    # println("rest: $(logdet(Σ) + d * log(2π))")
    try
        return val - (logdet(Σ) + d * log(2π)) / 2
    catch DomainError
        return 1e8
    end
end

approximate_evidence(x...) = exp(-approximate_negative_log_evidence(x...))

# undamped newton algorithm for fine-tuning
function newton!(f, θ::AbstractVector; min_step::Real = 1e-10, maxiter::Int = 5)
    N = OptimizationAlgorithms.SaddleFreeNewton(f, θ)
    d = OptimizationAlgorithms.direction(N, θ)
    i = 1
    while maximum(abs, d) > min_step * maximum(abs, θ) && i < maxiter
        OptimizationAlgorithms.update!(N, θ)
        d = N.d
        i += 1
    end
    return θ
end

# TODO: move these into optimize.jl in CrystalShift, to make everything cleaner
function ls_objective(PM::PhaseModel, log_θ::AbstractVector,
        x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
        std_noise::Real, mean_log_θ::AbstractVector, std_θ::AbstractVector)
    ls_residual(PM, log_θ, x, y, y_uncer, std_noise) + ls_regularizer(PM, log_θ, mean_log_θ, std_θ)
end

function ls_residual(PM::PhaseModel, log_θ::AbstractVector, x::AbstractVector,
        y::AbstractVector, y_uncer::AbstractVector, std_noise::Real)
    r = zeros(promote_type(eltype(log_θ), eltype(x), eltype(y)), length(x))
    r = _weighted_residual!(PM, log_θ, x, y, y_uncer, r, std_noise)
    return sum(abs2, r)
end

#################################### KL Divergence #############################
function kl_objective(phases::AbstractVector, log_θ::AbstractVector,
        x::AbstractVector, y::AbstractVector,
        mean_log_θ::AbstractVector, std_θ::AbstractVector, λ::Real)
    θ = exp.(log_θ)
    r_θ = evaluate!(phases, θ, x) # reconstruction of phases, IDEA: pre-allocate result (one for Dual, one for Float)
    r_θ ./= exp(1) # since we are not normalizing the inputs, this rescaling has the effect that kl(α*y, y) has the optimum at α = 1
    kl(r_θ, y) + λ * ls_regularizer(log_θ, mean_log_θ, std_θ)
end

function ls_regularizer(PM::PhaseModel, log_θ::AbstractVector, mean_log_θ::AbstractVector, std_θ::AbstractVector)
    p = zero(eltype(log_θ))
    bg_param_num = get_param_nums(PM.background)
    θ_cp = log_θ[1:end - bg_param_num ]
    θ_bg = log_θ[end - bg_param_num + 1 : end]
    @inbounds @simd for i in eachindex(θ_cp)
        p += ((log_θ[i] - mean_log_θ[i]) / (sqrt(2)*std_θ[i]))^2
    end
    p += _prior(PM.background, θ_bg)
    return p
end
