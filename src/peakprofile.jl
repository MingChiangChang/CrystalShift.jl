##################
# This code is taken from https://github.com/SebastianAment/PhaseMapping.jl
# With author's permission so remove the dependency of this module on PhaseMapping.jl
##################
abstract type PeakProfile{T} end

@inline locationscale(f, x, μ, σ) = f((x-μ)/σ)
(P::PeakProfile)(x::Real, μ::Real, σ::Real) = locationscale(P, x, μ, σ)
(P::PeakProfile)(x::Real, c::Real, μ::Real, σ::Real) = c*P(x, μ, σ)

get_param_nums(P::PeakProfile) = 0
get_free_params(P::PeakProfile) = []

############################### Lorentzian function ##############################
struct LorentzianProfile{T} <: PeakProfile{T} end
const Lorentz = LorentzianProfile
Lorentz() = Lorentz{Float64}()
(P::Lorentz)(x::Real) = inv(1 + x^2)

############################### Gaussian function ##############################
struct GaussianProfile{T} <: PeakProfile{T} end
const Gauss = GaussianProfile
Gauss() = Gauss{Float64}()
function Gauss(a::AbstractArray) 
    isempty(a) && error("Gauss recieves a parameter")
    Gauss()
end
(P::Gauss)(x::Real) = special_exp(-x^2/2)

@inline sigmoid(a::Real) = 1/(1+exp(-a))
@inline inverse_sig(a::Real) = -log((1-a)/a)

function inverse_sig()
    return []
end

############################# pseudo-voigt function ############################
struct PseudoVoigtProfile{T} <: PeakProfile{T}
   α::T
   sig_α::T

   function PseudoVoigtProfile{T}(a::T) where T<:Real
      new{T}(a, sigmoid(a))
   end
end

const PseudoVoigt = PseudoVoigtProfile
PseudoVoigt(a::Float64) = PseudoVoigt{Float64}(a)
(P::PseudoVoigt)(x::Real) = (-0.5+P.sig_α) * Lorentz()(x) + (1.5-P.sig_α) * Gauss()(x)
get_param_nums(P::PseudoVoigtProfile) = length(P.α)
get_free_params(P::PseudoVoigtProfile) = [P.sig_α] # This is in sigmoid space!!

function PseudoVoigt(a::AbstractVector) 
   length(a) == 1 || error("PseudoVoigt receive more than one params")
   PseudoVoigt{eltype(a)}(a[1])
end

struct FixedPseudoVoigtProfile{T} <: PeakProfile{T}
   α::T
end

const FixedPseudoVoigt = FixedPseudoVoigtProfile
FixedPseudoVoigt(a::Float64) = FixedPseudoVoigt{Float64}(a)
(P::FixedPseudoVoigt)(x::Real) = P.α * Lorentz()(x) + (1-P.α) * Gauss()(x)
get_param_nums(P::FixedPseudoVoigt) = 0
get_free_params(P::FixedPseudoVoigt) = []

########################## mixture of peak profiles ############################
# θ parameters of mixture (peak parameters x number of peaks)
# function mixture!(y::AbstractVector, P::PeakProfile,
#                   x::AbstractVector, θ::AbstractMatrix)
#    for j in 1:size(θ, 2)
#       @simd for i in eachindex(x)
#          @inbounds y[i] = P(x[i], θ[:, j]...)
#       end
#    end
#    return y
# end

# function mixture(P::PeakProfile, x::AbstractVector, θ::AbstractMatrix)
#    mixture!(similar(x), P, x, θ)
# end

# # θ parameters of mixture (peak parameters x number of peaks x number of patterns)
# function mixture!(Y::AbstractMatrix, P::PeakProfile,
#                   x::AbstractVector, θ::AbstractArray{<:Real, 3})
#    k = size(θ, 3)
#    size(Y, 1) == length(x) || throw(DimensionMismatch(""))
#    size(Y, 2) == k || throw(DimensionMismatch(""))
#    for i in 1:k
#       y_i = @view Y[:, i]
#       θ_i = @view θ[:, :, i]
#       mixture!(y_i, P, x, θ_i)
#    end
#    return Y
# end

# function mixture(P::PeakProfile, x::AbstractVector, θ::AbstractArray{<:Real, 3})
#    Y = zeros(eltype(x), (length(x), size(θ, 3)))
#    mixture!(Y, P, x, θ)
# end


# ###### gradients
# # .031ns vs 11.404 ns (lorentz vs gauss)
# @inline function gradient(::Lorentz, x::Real, c::Real, μ::Real, σ::Real)
#    ξ = (x-μ)/σ
#    e = inv(1 + ξ^2)
#    y = c * e # value
#    dc = e # this stays the same regardless of peak function
#    dμ = y * e * 2ξ/σ
#    dσ = dμ * ξ
#    return y, (dc, dμ, dσ)
# end

# # gauss gradient definition
# @inline function gradient(::Gauss, x::Real, c::Real, μ::Real, σ::Real)
#    ξ = (x-μ)/σ
#    e = exp(-ξ^2/2)
#    y = c * e # value
#    dc = e
#    yξ = y*ξ/σ # temporary to pool computation
#    dμ = yξ
#    dσ = yξ*ξ
#    return y, (dc, dμ, dσ)
# end
