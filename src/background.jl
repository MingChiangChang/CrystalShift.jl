using CovarianceFunctions
using LinearAlgebra
const DEFAULT_RANK_TOL = 1e-6

struct BackgroundModel{T, KT, AT<:AbstractMatrix{T}, UT<:AbstractMatrix, ST<:AbstractVector, LT, CT}
    k::KT # kernel function
    K::AT # kernel matrix
    U::UT # low rank approximation
    S::ST # singular values of SVD used to construct low rank approximation to K
    λ::LT # coefficient for prior over model parameters c
    c::CT # temporary storage for model parameters
end

get_param_nums(B::BackgroundModel) = length(B.S)
get_free_params(B::BackgroundModel) = B.c

function BackgroundModel(x::AbstractVector, k, l::Real, λ::Real = 1; rank_tol::Real = DEFAULT_RANK_TOL)
    k = CovarianceFunctions.Lengthscale(k, l)
    K, U, S = low_rank_background_model(x, k, rank_tol)
    c = zeros(length(S))
    BackgroundModel(k, K, U, S, λ, c)
end

function reconstruct_BG!(θ::AbstractVector, B::BackgroundModel)
    param_num = get_param_nums(B)
    new_B = BackgroundModel(B.k, B.K, B.U, B.S, B.λ, θ[1:param_num])
    return θ[param_num+1:end], new_B
end
# function BackgroundModel(BG::BackgroundModel,  c::)

reconstruct_BG!(θ::AbstractVector, B::Nothing) = θ, nothing

LinearAlgebra.svd(G::Gramian) = svd(Matrix(G))

# x is vector of q values in case of XRD
# k is kernel function
# rank_tol is the singular value threshold below which we truncate and SVD to approximate the full kernel matrix
# Output: Ui is the low rank matrix given coefficients α, Ui*α is the background model.
# To be true to the original kernel model, have to add least-squares regularizer on α
# with variances equal to the singular values Si.
function low_rank_background_model(x::AbstractVector, k, rank_tol::Real = DEFAULT_RANK_TOL)
    K = gramian(k, x)
    U, S, V = svd(K)
    i = findfirst(<(rank_tol), S)
    Ui = @view U[:, 1:i]
    Si = @view S[1:i]
    return K, Ui, Si
end

# function evaluate(B::BackgroundModel, c::AbstractVector)
#     y = zeros(eltype(c), size(B.K, 1))
#     evaluate!(y, B, c)
# end


# function evaluate!(y::AbstractVector, B::BackgroundModel, c::AbstractVector)
#     A = isnothing(B.U) ? B.K : B.U
#     mul!(y, A, c)
# end

function evaluate!(y::AbstractVector, B::BackgroundModel, θ::AbstractVector, x::AbstractVector)
    reconstruct_BG!(θ, B)
    evaluate!(y, B, x)
end

function evaluate(B::BackgroundModel, x::AbstractVector)
    y = zero(x)
    evaluate!(y, B, x)
end

function evaluate(B::BackgroundModel, θ::AbstractVector, x::AbstractVector)
    y = zero(x)
    evaluate!(y, B, θ, x)
end

function evaluate!(y::AbstractVector, B::BackgroundModel, x::AbstractVector)
    A = isnothing(B.U) ? B.K : B.U
    y .+= mul!(zero(x), A, B.c)
    y
end

function evaluate_residual!(BG::BackgroundModel, θ::AbstractVector,
                            x::AbstractVector, r::AbstractVector)
    evaluate_residual!(reconstruct_BG(BG, θ), x, r)
end

function evaluate_residual!(B::BackgroundModel, x::AbstractVector, r::AbstractVector)
    A = isnothing(B.U) ? B.K : B.U
    y = zero(r)
    r .-= mul!(y, A, B.c)
    r
end

function _prior(B::BackgroundModel, c::AbstractVector)
    p = zero(c)
    lm_prior!(p, B, c)
    sum(abs2, p)
end

_prior(B::Nothing, c::AbstractVector) = 0

# least-squares prior for background model to be used in conjunction with LevenbergMarquart
# computes
function lm_prior!(p::AbstractVector, B::BackgroundModel, c::AbstractVector)
    if isnothing(B.U)
        @. p = B.λ * c
    else
        @. p = B.λ * c / B.S
    end
end

function lm_prior!(p::AbstractVector, B::BackgroundModel)
    if isnothing(B.U)
        @. p = B.λ * B.c
    else
        @. p = B.λ * B.c / B.S
    end
end

function lm_prior!(p::AbstractVector, B::Nothing, c::AbstractVector) end
