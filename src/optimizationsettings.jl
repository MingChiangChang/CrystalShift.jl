const ALLOWED_OBJECTIVE = ["LS", "KL"]
const PRIOR_LENGTH = 3
const DEFAULT_TOL = 1e-8
@exported_enum OptimizationMethods LM Newton

function extend_priors(mean_θ::AbstractVector, std_θ::AbstractVector,
    phases::AbstractVector{<:CrystalPhase})
    extend_priors(mean_θ, std_θ, [phase.cl for phase in phases])
end

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
            full_std_θ[start + n + 2] = Inf
            start += (n+3)
        else
            start += (n+2)
        end
    end
    return full_mean_θ, full_std_θ
end

struct Priors{T}
    std_noise::Real
    mean_θ::AbstractVector{T}
    std_θ::AbstractVector{T}

    function Priors{V}(phases::AbstractVector{<:CrystalPhase}, std_noise::Real,
                       mean_θ::AbstractVector{V}, std_θ::AbstractVector{V}) where V<:Real
        length(mean_θ) == length(std_θ) == PRIOR_LENGTH || error("Prior must have length of $(PRIOR_LENGTH)")
        std_noise > 0 || error("std_noise must be larger than zero")
        
        mean_θ, std_θ = extend_priors(mean_θ, std_θ, phases)
        new{V}(std_noise, mean_θ, std_θ)
    end
end


Priors{T}() where T<:Real = Priors{T}(0.1, [1., 1e-4, .2], [.02, 100, 1.])

struct OptimizationSettings{T}
    priors::Priors{T}
    maxiter::Int
    regularization::Bool
    method::OptimizationMethods
    objective::String
    verbose::Bool
    tol::Float64

    function OptimizationSettings{V}(priors::Priors{V}, maxiter::Int, regularization::Bool,
                                     method::OptimizationMethods, objective::String,
                                     verbose::Bool, tol::Float64) where V<:Real
        maxiter > 0 || error("max_iter must be > 0")
        objective in ALLOWED_OBJECTIVE || ("Objective string not in allowed objective")
        new{V}(priors, maxiter, regularization, method, objective, verbose, tol)
    end
end

function OptimizationSettings{V}(phases::AbstractVector{<:CrystalPhase},
                                 std_noise::Real, mean_θ::AbstractVector{V}, std_θ::AbstractVector{V}, 
                                 maxiter::Int, regularization::Bool, method::OptimizationMethods, 
                                 objective::String, verbose::Bool, tol::Float64) where V<:Real
    pr = Priors{V}(phases, std_noise, mean_θ, std_θ)
    OptimizationSettings{V}(pr, maxiter, regularization, method, objective, verbose, tol)
end