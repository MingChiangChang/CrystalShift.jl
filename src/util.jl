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
