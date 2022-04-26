struct Wildcard{K, VK<:AbstractVector{K}, PR, T} <: AbstractPhase
    peak_locs::VK
    activations::VK
    name::String
    σ::VK
    profile::PR

    mean_θ::T
    std_θ::T
end

function Wildcard(peak_locs::AbstractVector, activations::AbstractVector,
                  name::String, σ::AbstractVector,
                  profile::PeakProfile, std_θ::AbstractVector)
    mean_θ = log.([peak_locs..., activations..., σ...])
    Wildcard(peak_locs, activations, name, σ, profile, mean_θ, std_θ)
end

function Wildcard(W::Wildcard, θ::AbstractVector)
    l = length(W.peak_locs)
    Wildcard(θ[1:l], θ[l+1:2*l], W.name, θ[2*l+1:3*l], W.profile, W.mean_θ, W.std_θ)
end

get_param_nums(W::Wildcard) = 3*length(W.peak_locs)
get_free_params(W::Wildcard) = [W.peak_locs..., W.activations..., W.σ...]

function evaluate!(y::AbstractVector, W::Wildcard, x::AbstractVector)
    @simd for i in eachindex(W.peak_locs)
        @. y +=  W.activations[i] * W.profile((x - W.peak_locs[i]) / W.σ[i])
    end
    y
end

function evaluate(W::Wildcard, x::AbstractVector)
    y = zero(x)
    evaluate!(y, W, x)
end

function evaluate_residual!(W::Wildcard, x::AbstractVector, r::AbstractVector)
    @simd for i in eachindex(W.peak_locs)
        @. r -=  W.activations[i] * W.profile((x - W.peak_locs[i]) / W.σ[i])
    end
    r
end

function evaluate!(y::AbstractVector, W::Wildcard, θ::AbstractVector,
                   x::AbstractVector)
    Wildcard(W, θ)(x, y)
end

function lm_prior!(p::AbstractVector, W::AbstractVector{<:Wildcard},
                   θ::AbstractVector)
    start = 1
    for i in eachindex(W)
        p .= (θ[start : start+get_param_nums(W[i])-1] .- W[i].mean_θ) ./ (sqrt(2) .* W[i].std_θ)
        start += get_param_nums(W[i])
    end
    p
end