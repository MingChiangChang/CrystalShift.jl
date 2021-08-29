cast(f::AbstractVector, type::Type) = map(x->parse(type, x), f)

# "chess_x_y_tpeak_dwell"
function parse_cond(cond::AbstractString, T::Type)
    conds = split(cond, "_")[2:end]
    return cast(conds, T)
end
