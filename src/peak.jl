struct Peak{T<:Real, H<:Real}
    h::T # Stored separately to be explicit
    k::T
    l::T

    I::H # Intensity
end

# Fallback function and for triclinic
function (cl::Crystal)(P::Peak)
    2pi/volume(cl)*sqrt(P.h^2)
end
