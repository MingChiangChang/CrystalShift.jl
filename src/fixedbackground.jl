struct FixedBackground{T, WT<:AbstractVector{T}, V, W, Z} <: AbstractBackground
    basis::WT
    a::V # Linear scaling factor
    c::W # Constant background factor
    λ::Z # Regularization factor
end

get_free_params(FBG::FixedBackground) = [FBG.a, FBG.c]
get_param_nums(FBG::FixedBackground) = 2


function reconstruct_BG!(θ::AbstractVector, B::FixedBackground)
    new_B = FixedBackground(B.basis, θ[1], θ[2], B.λ)
    return θ[3:end], new_B
end

function evaluate!(y::AbstractVector, B::FixedBackground, x::AbstractVector)
    @. y += B.a * B.basis + B.c # This will fail if basis does not have the right dimension
    y
end

function evaluate_residual!(B::FixedBackground, x::AbstractVector, r::AbstractVector)
    @. r -= B.a * B.basis + B.c
    r
end

function _prior(B::FixedBackground, c::AbstractVector)
    p = zero(c)
    lm_prior!(p, B, c)
    sum(abs2, p)
end

function lm_prior!(p::AbstractVector, B::FixedBackground, c::AbstractVector)
    @. p = B.λ * c
end

function lm_prior!(p::AbstractVector, B::FixedBackground)
    @. p = B.λ * get_free_params
end
