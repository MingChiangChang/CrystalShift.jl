cast(f::AbstractVector, type::Type) = map(x->parse(type, x), f)

# "chess_x_y_tpeak_dwell"
function parse_cond(cond::AbstractString, T::Type)
    conds = split(cond, "_")[2:end]
    return cast(conds, T)
end

# TODO Should have a better way to identify center
function get_weighted_center(H::AbstractMatrix)
    avg = []
    for act in eachrow(H)
        ind = collect(1:size(act)[1])
        s = sum(act)
        if s>0
            push!(avg, Int(floor( sum(act.*ind)/sum(act) )))
        else
            print(size(act)[1]/2)
            push!(avg, Int(floor(size(act)[1]/2)))
        end
    end
    avg
end

## convenient macros
macro exported_enum(name, args...)
    esc(quote
        @enum($name, $(args...))
        export $name
        $([:(export $arg) for arg in args]...)
        end)
end

########################## objective functions #################################
# Kullback-Leibler divergence
function kl(p::Real, q::Real)
    !(isnan(p) || isinf(p) || isnan(q) || isinf(q)) || throw(DomainError("Nan or inf input"))
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
