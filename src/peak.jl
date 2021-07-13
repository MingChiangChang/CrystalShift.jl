struct Peak{T<:Real, H<:Real}
    h::T # Stored separately to be explicit
    k::T
    l::T

    I::H # Intensity
end
