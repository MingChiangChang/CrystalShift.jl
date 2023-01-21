const ALLOWED_OBJECTIVE = ["LS", "KL"]
const PEAK_PRIOR_LENGTH = 1
const DEFAULT_TOL = 1e-8
@exported_enum OptimizationMethods LM Newton bfgs l_bfgs
@exported_enum OptimizationMode Simple EM WithUncer

# function extend_priors(mean_θ::AbstractVector, std_θ::AbstractVector,
#     phases::AbstractVector{<:AbstractPhase})
#     extend_priors(mean_θ, std_θ, [phase.cl for phase in phases])
# end

struct Priors{T}
    std_noise::Real
    mean_θ::AbstractVector{T}
    std_θ::AbstractVector{T}

    function Priors{V}(std_noise::Real, mean_θ::AbstractVector{V}, std_θ::AbstractVector{V}) where V<:Real
        length(mean_θ) == length(std_θ) || error("Prior must have same length")
        std_noise > 0 || error("std_noise must be larger than zero")

        new{V}(std_noise, mean_θ, std_θ)
    end

    # function Priors{V}(phases::Nothing, std_noise::Real,
    #     mean_θ::AbstractVector{V}, std_θ::AbstractVector{V}) where V<:Real
    #     new{V}(std_noise, [], [])
    # end
end

extend_priors(pr::Priors, phases::AbstractVector{<:AbstractPhase}) = extend_priors(pr.mean_θ, pr.std_θ, phases)
extend_priors(pr::Priors, phase_model::PhaseModel) = extend_priors(pr, phase_model.CPs)
extend_priors(pr::Priors, phases::Nothing) = [], []

function extend_priors(mean_θ::AbstractVector, std_θ::AbstractVector,
                        phases::AbstractVector{<:CrystalPhase})
    totl_params = sum([get_param_nums(phase) for phase in phases])
    full_mean_θ = zeros(totl_params)
    full_std_θ = zeros(totl_params)
    start = 1
    for phase in phases
        n = phase.cl.free_param
        full_mean_θ[start:start+n-1] = mean_θ[1].*get_free_lattice_params(phase)
        full_std_θ[start:start+n-1] = std_θ[1].*get_free_lattice_params(phase)#repeat(std_θ[1, :], n)
        full_mean_θ[start + n:start + n + 1] = mean_θ[2:3]
        full_std_θ[start + n:start + n + 1] = std_θ[2:3]
        if phase.profile isa PseudoVoigt
            full_mean_θ[start + n + 2] = 0.5
            full_std_θ[start + n + 2] = 10.
        end
        start += get_param_nums(phase)
    end
    return full_mean_θ, full_std_θ
end

function extend_priors(mean_θ::AbstractVector, std_θ::AbstractVector,
                       phases::Nothing)
    return [], []
end

function extend_priors(mean_θ::AbstractVector, std_θ::AbstractVector,
        phases::AbstractVector{<:PeakModCP})
    totl_params = sum([get_param_nums(phase) for phase in phases])
    full_mean_θ = zeros(totl_params)
    full_std_θ = zeros(totl_params)
    start = 1

    for phase in phases
        n = get_param_nums(phase)
        full_mean_θ[start:start+n-1] = get_free_params(phase) .* mean_θ
        full_std_θ[start:start+n-1] = repeat(std_θ[1:1], n)
        start += n
    end

    full_mean_θ, full_std_θ
end




Priors{T}() where T<:Real = Priors{T}(0.1, [1., 1e-4, .2], [.02, 100, 1.])

struct OptimizationSettings{T}
    priors::Priors{T}
    maxiter::Int
    regularization::Bool
    method::OptimizationMethods
    objective::String
    optimize_mode::OptimizationMode
    em_loop_num::Integer
    verbose::Bool
    tol::Float64

    function OptimizationSettings{V}(priors::Priors{V},
                                     maxiter::Int = 128,
                                     regularization::Bool =true,
                                     method::OptimizationMethods = LM,
                                     objective::String = "LS",
                                     optimize_mode::OptimizationMode = Simple,
                                     em_loop_num::Integer=8,
                                     verbose::Bool=false,
                                     tol::Float64=DEFAULT_TOL) where V<:Real
        maxiter > 0 || error("max_iter must be > 0")
        objective in ALLOWED_OBJECTIVE || ("Objective string not in allowed objective")
        new{V}(priors, maxiter, regularization, method, objective, optimize_mode, em_loop_num, verbose, tol)
    end
end

function OptimizationSettings{V}(
                                 std_noise::Real, mean_θ::AbstractVector{V}, std_θ::AbstractVector{V},
                                 maxiter::Int = 128,
                                 regularization::Bool =true,
                                 method::OptimizationMethods = LM,
                                 objective::String = "LS",
                                 optimize_mode::OptimizationMode = Simple,
                                 em_loop_num::Integer=8,
                                 verbose::Bool=false,
                                 tol::Float64=DEFAULT_TOL) where V<:Real
    pr = Priors{V}(std_noise, mean_θ, std_θ)
    OptimizationSettings{V}(pr, maxiter, regularization, method, objective, optimize_mode, em_loop_num, verbose, tol)
end

function OptimizationSettings{Float64}()
    pr = Priors{Float64}()
    OptimizationSettings{Float64}(pr, 128, true, LM, "LS", Simple, 8, false, DEFAULT_TOL)
end
