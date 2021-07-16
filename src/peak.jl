struct Peak{T<:Real, H<:Real}
    h::T # Stored separately to be explicit
    k::T
    l::T

    I::H # Intensity
end

function (cl::Crystal)(P::Peak)
    2pi/volume(cl)*sqrt(cl.h^2)
end
