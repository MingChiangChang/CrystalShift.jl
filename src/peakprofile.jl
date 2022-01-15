abstract type PeakProfile{T} end

@inline locationscale(f, x, μ, σ) = f((x-μ)/σ)
(P::PeakProfile)(x::Real, μ::Real, σ::Real) = locationscale(P, x, μ, σ)
(P::PeakProfile)(x::Real, c::Real, μ::Real, σ::Real) = c*P(x, μ, σ)

############################### Lorentzian function ##############################
struct LorentzianProfile{T} <: PeakProfile{T} end
const Lorentz = LorentzianProfile
Lorentz() = Lorentz{Float64}()
(P::Lorentz)(x::Real) = inv(1 + x^2)

############################### Gaussian function ##############################
struct GaussianProfile{T} <: PeakProfile{T} end
const Gauss = GaussianProfile
Gauss() = Gauss{Float64}()
(P::Gauss)(x::Real) = exp(-x^2/2)

############################# pseudo-voigt function ############################
struct PseudoVoigtProfile{T} <: PeakProfile{T}
   α::T
end
const PseudoVoigt = PseudoVoigtProfile
(P::PseudoVoigt)(x::Real) = P.α * Lorentz()(x) + (1-P.α) * Gauss()(x)

########################## mixture of peak profiles ############################
# θ parameters of mixture (peak parameters x number of peaks)
function mixture!(y::AbstractVector, P::PeakProfile,
                  x::AbstractVector, θ::AbstractMatrix)
   for j in 1:size(θ, 2)
      @simd for i in eachindex(x)
         @inbounds y[i] = P(x[i], θ[:, j]...)
      end
   end
   return y
end

function mixture(P::PeakProfile, x::AbstractVector, θ::AbstractMatrix)
   mixture!(similar(x), P, x, θ)
end

# θ parameters of mixture (peak parameters x number of peaks x number of patterns)
function mixture!(Y::AbstractMatrix, P::PeakProfile,
                  x::AbstractVector, θ::AbstractArray{<:Real, 3})
   k = size(θ, 3)
   size(Y, 1) == length(x) || throw(DimensionMismatch(""))
   size(Y, 2) == k || throw(DimensionMismatch(""))
   for i in 1:k
      y_i = @view Y[:, i]
      θ_i = @view θ[:, :, i]
      mixture!(y_i, P, x, θ_i)
   end
   return Y
end

function mixture(P::PeakProfile, x::AbstractVector, θ::AbstractArray{<:Real, 3})
   Y = zeros(eltype(x), (length(x), size(θ, 3)))
   mixture!(Y, P, x, θ)
end


###### gradients
# .031ns vs 11.404 ns (lorentz vs gauss)
@inline function gradient(::Lorentz, x::Real, c::Real, μ::Real, σ::Real)
   ξ = (x-μ)/σ
   e = inv(1 + ξ^2)
   y = c * e # value
   dc = e # this stays the same regardless of peak function
   dμ = y * e * 2ξ/σ
   dσ = dμ * ξ
   return y, (dc, dμ, dσ)
end

# gauss gradient definition
@inline function gradient(::Gauss, x::Real, c::Real, μ::Real, σ::Real)
   ξ = (x-μ)/σ
   e = exp(-ξ^2/2)
   y = c * e # value
   dc = e
   yξ = y*ξ/σ # temporary to pool computation
   dμ = yξ
   dσ = yξ*ξ
   return y, (dc, dμ, dσ)
end

# (P::PseudoVoigt)(x::Real, α::Real) = α * Lorentz()(x) + (1-α) * Gauss()(x)
# (P::PseudoVoigt)(x::Real, μ::Real, σ::Real, α::Real) = locationscale(P, x, μ, σ, α)
# (P::PseudoVoigt)(x::Real, c::Real, μ::Real, σ::Real, α::Real) = c*P(x, μ, σ, α)

# @inline function gradient(::PseudoVoigt, x::Real, c::Real, μ::Real, σ::Real, α::Real)
#    yl, ∇l = gradient(Lorentz(), x, c, μ, σ)
#    yg, ∇g = gradient(Gauss(), x, c, μ, σ)
#    y = α * yl + (α-1) * yg
#    ∇ = (begin @. α * ∇l + (α-1) * ∇g end, (yl+yg)/2)
#    return y, ∇
# end
