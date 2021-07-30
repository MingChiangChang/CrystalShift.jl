struct Peak{T<:Real, H<:Real}
    h::T # Stored separately to be explicit
    k::T
    l::T

    q::T # Reference for now, probably don't need in final version
    I::H # Intensity
end

function Peak(s::String)
    info = split(s, ',')
    println(info)
    length(info) == 5 || throw("info must has length of 5")
    h, k, l, q, I = cast(info, Float64)
    Peak(h, k, l, q, I)
end
