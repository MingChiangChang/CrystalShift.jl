cast(f::AbstractVector, type::Type) = map(x->parse(type, x), f)

# "chess_x_y_tpeak_dwell"
function parse_cond(cond::AbstractString, T::Type)
    conds = split(cond, "_")[4:end]
    return cast(conds, T)
end

# TODO Should have a better way to identify center
function get_weighted_center(H::AbstractMatrix)
    avg = []
    for act in eachrow(H)
        ind = collect(1:size(act, 1))
        s = sum(act)
        if s>0
            push!(avg, Int(floor( sum(act.*ind)/sum(act) )))
        else
            print(size(act)[1]/2)
            push!(avg, Int(floor(size(act)[1]/2)))
        end
    end
    sum(avg)/length(avg)
end

## convenient macros
macro exported_enum(name, args...)
    esc(quote
        @enum($name, $(args...))
        export $name
        $([:(export $arg) for arg in args]...)
        end)
end

function volume(a::Real, b::Real, c::Real,
                α::Real, β::Real, γ::Real)
        return (a * b * c
        * sqrt( 1+2*cos(α)*cos(β)*cos(γ)
              - cos(α)^2 - cos(β)^2 - cos(γ)^2 ) )
end

function special_inv(x::Real)
    if x > 100.
        return 1E-12
    else
        return inv(x)
    end
end

function special_exp(x::Real)
    if x < -10.
        return 1E-12
    else
        return exp(x)
    end
end



########################## objective functions #################################
# Kullback-Leibler divergence
function kl(p::Real, q::Real)
    !(isnan(p) || isinf(p) || isnan(q) || isinf(q)) || throw(DomainError("$(p) Nan or inf input"))
    (p > 0 && q > 0) || throw(DomainError("kl divergence undefined for negative inputs"))
    iszero(p) ? zero(p) : (p*log(p) - p*log(q)) # -p*log(q/p) # p*log(p/q) # is equal to
end

function kl(P::AbstractArray, Q::AbstractArray)
    length(P) == length(Q) || throw(DimensionMismatch("Lengths of inputs arrays do not match"))
    val = zero(promote_type(eltype(P), eltype(Q)))
    @inbounds @simd for i in eachindex(P)
        val += kl(P[i], Q[i])
    end
    return val
end


# negative Poisson likelihood, generalized for real (not just integer k)
function negative_log_poisson(λ::Real, k::Real)
    λ + loggamma(k) - k*log(λ)
end

# if we are optimizating w.r.t. λ, gamma(k) is irrelevant, so we remove it for this function
# NOTE: this is equivalent to optimizing kl w.r.t. 2nd argument and having l1 penalty
function negative_log_poisson_λ(λ::Real, k::Real)
    λ - k*log(λ)
end

poisson(λ::Real, k::Real) = λ^k * exp(-λ) / gamma(k)
