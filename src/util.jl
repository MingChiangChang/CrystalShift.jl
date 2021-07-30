cast(f::AbstractVector, type::Type) = map(x->parse(type, x), f)
