# Common Rietveld background models that are linear in their coefficients:
# y_bg = basis * c, where basis is precomputed on the grid x at construction.
# Coefficients are refined jointly with the phases (not log-transformed),
# with an independent least-squares prior λ[i] * c[i] per coefficient.
struct BasisBackground{T, MT<:AbstractMatrix{T}, CT<:AbstractVector, LT<:AbstractVector} <: AbstractBackground
    basis::MT # length(x) × number of coefficients
    c::CT     # coefficients
    λ::LT     # per-coefficient regularization, 0 means unregularized
end

function BasisBackground(basis::AbstractMatrix, λ::Union{Real, AbstractVector} = 0.)
    n = size(basis, 2)
    λ = λ isa Real ? fill(float(λ), n) : collect(float.(λ))
    length(λ) == n || throw(ArgumentError("length(λ) = $(length(λ)) does not match number of basis functions $n"))
    BasisBackground(basis, zeros(n), λ)
end

is_linear(::BasisBackground) = true
background_basis(B::BasisBackground, x::AbstractVector) = B.basis
get_param_nums(B::BasisBackground) = size(B.basis, 2)
get_free_params(B::BasisBackground) = B.c

function reconstruct_BG!(θ::AbstractVector, B::BasisBackground)
    n = get_param_nums(B)
    return θ[n+1:end], BasisBackground(B.basis, θ[1:n], B.λ)
end

function evaluate!(y::AbstractVector, B::BasisBackground, x::AbstractVector)
    check_basis_size(B, y)
    y .+= B.basis * B.c
    y
end

function evaluate_residual!(B::BasisBackground, x::AbstractVector, r::AbstractVector)
    check_basis_size(B, r)
    r .-= B.basis * B.c
    r
end

function check_basis_size(B::BasisBackground, y::AbstractVector)
    size(B.basis, 1) == length(y) || throw(DimensionMismatch("background basis was built for $(size(B.basis, 1)) points but got $(length(y))"))
end

function _prior(B::BasisBackground, c::AbstractVector)
    p = zero(c)
    lm_prior!(p, B, c)
    sum(abs2, p)
end

lm_prior!(p::AbstractVector, B::BasisBackground, c::AbstractVector) = (@. p = B.λ * c)
lm_prior!(p::AbstractVector, B::BasisBackground) = (@. p = B.λ * B.c)

# maps x linearly onto [-1, 1]
function rescale_to_unit(x::AbstractVector)
    lo, hi = extrema(x)
    hi > lo || throw(ArgumentError("x must span a nonzero range"))
    @. 2 * (x - lo) / (hi - lo) - 1
end

# λ for the constant term is λ0 (default unregularized), λ for the rest is λ
coefficient_priors(n::Int, λ::Real, λ0::Real) = [float(λ0); fill(float(λ), n - 1)]

"""
    PolynomialBackground(x, order; λ = 1., λ0 = 0., inverse = false)

Power series Σ cᵢ tⁱ for i = 0..order, where t is x rescaled to [-1, 1]
for conditioning. `inverse = true` appends a term min(x)/x, which models
low-angle air scattering. Prefer `ChebyshevBackground` for order ≳ 6.
"""
function PolynomialBackground(x::AbstractVector, order::Int; λ::Real = 1., λ0::Real = 0., inverse::Bool = false)
    order >= 0 || throw(ArgumentError("order must be non-negative"))
    t = rescale_to_unit(x)
    basis = reduce(hcat, [t .^ i for i in 0:order])
    λs = coefficient_priors(order + 1, λ, λ0)
    if inverse
        all(>(0), x) || throw(ArgumentError("inverse term requires x > 0"))
        basis = hcat(basis, minimum(x) ./ x)
        push!(λs, float(λ))
    end
    BasisBackground(basis, λs)
end

"""
    ChebyshevBackground(x, order; λ = 1., λ0 = 0.)

Chebyshev polynomials of the first kind Σ cᵢ Tᵢ(t) for i = 0..order,
where t is x rescaled to [-1, 1] (GSAS background type 1).
"""
function ChebyshevBackground(x::AbstractVector, order::Int; λ::Real = 1., λ0::Real = 0.)
    order >= 0 || throw(ArgumentError("order must be non-negative"))
    t = rescale_to_unit(x)
    basis = ones(length(x), order + 1)
    order >= 1 && (basis[:, 2] .= t)
    for i in 3:order+1
        @. basis[:, i] = 2t * basis[:, i-1] - basis[:, i-2]
    end
    BasisBackground(basis, coefficient_priors(order + 1, λ, λ0))
end

"""
    CosineBackground(x, order; λ = 1., λ0 = 0.)

Cosine Fourier series Σ cᵢ cos(iπu) for i = 0..order, where u is x rescaled
to [0, 1] (in the spirit of GSAS background type 2).
"""
function CosineBackground(x::AbstractVector, order::Int; λ::Real = 1., λ0::Real = 0.)
    order >= 0 || throw(ArgumentError("order must be non-negative"))
    u = (rescale_to_unit(x) .+ 1) ./ 2
    basis = reduce(hcat, [cos.(i * π .* u) for i in 0:order])
    BasisBackground(basis, coefficient_priors(order + 1, λ, λ0))
end

"""
    InterpolationBackground(x, anchors; λ = 0.)
    InterpolationBackground(x, n::Int; λ = 0.)

Piecewise-linear background through refinable heights at the `anchors`
(GSAS type 7 / FullProf background points). Constant beyond the outermost
anchors. The integer form places `n` anchors evenly across the range of x.
"""
function InterpolationBackground(x::AbstractVector, anchors::AbstractVector; λ::Real = 0.)
    a = sort(collect(float.(anchors)))
    length(a) >= 2 || throw(ArgumentError("need at least two anchors"))
    allunique(a) || throw(ArgumentError("anchors must be distinct"))
    basis = zeros(length(x), length(a))
    for (k, xk) in enumerate(x)
        if xk <= a[1]
            basis[k, 1] = 1
        elseif xk >= a[end]
            basis[k, end] = 1
        else
            j = searchsortedlast(a, xk)
            w = (xk - a[j]) / (a[j+1] - a[j])
            basis[k, j] = 1 - w
            basis[k, j+1] = w
        end
    end
    BasisBackground(basis, λ)
end

function InterpolationBackground(x::AbstractVector, n::Int; λ::Real = 0.)
    InterpolationBackground(x, range(extrema(x)..., length = n); λ = λ)
end
