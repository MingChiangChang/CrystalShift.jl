using PhaseMapping
using DelimitedFiles

# Abstract type to as supertype for 7 different crystal systems
abstract type Crystal end

# TODO use Base.getproperty to make things easier
struct Triclinic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param = 6

    function Triclinic{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:AbstractFloat}
        new{T}(a, b, c, α, β, γ)
    end
end

struct Monoclinic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param = 3

    function Monoclinic{T}(a::T, c::T, β::T) where {T<:AbstractFloat}
        new{T}(a, a, c, pi/2, β, pi/2)
    end
end

struct Orthorhombic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param = 3
    function Orthorhombic{T}(a::T, b::T, c::T) where {T<:AbstractFloat}
        new{T}(a, b, c, pi/2, pi/2, pi/2)
    end
end

struct Tetragonal{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param = 2

    function Tetragonal{T}(a::T, c::T) where {T<:AbstractFloat}
        new{T}(a, a, c, pi/2, pi/2, pi/2)
    end
end

struct Rhombohedral{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param = 1

    function Rhombohedral{T}(a) where {T<:AbstractFloat}
        new{T}(a, a, a, pi/2, pi/2, 2*pi/3)
    end
end

struct Hexagonal{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param = 2

    function Hexagonal{T}(a::T, c::T) where {T<:AbstractFloat}
        new{T}(a, a, c, pi/2, pi/2, 2*pi/3)
    end
end

struct Cubic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param = 1

    # get_property to get b,c, α, β, γ
    function Cubic{T}(a::T) where {T<:AbstractFloat}
        new{T}(a, a, a, pi/2, pi/2, pi/2)
    end
end

check_not_equal(x...) = length(Set(x)) == length(x)
check_equal(x...) = all(y->y==x[1], x)

function isCubic(a::Real, b::Real, c::Real,
                 α::Real, β::Real, γ::Real)
    return check_equal(a, b, c) && check_equal(α, β, γ, pi/2)
end

function isTetragonal(a::Real, b::Real, c::Real,
                      α::Real, β::Real, γ::Real)
    return check_not_equal(a, c) && check_equal(a, b) &&\
           check_equal(a, b) && check_equal(α, β, γ, pi/2)
end

function isHexagonal(a::Real, b::Real, c::Real,
                      α::Real, β::Real, γ::Real)
    return check_not_equal(a, c) && check_equal(a, b) &&\
           check_equal(α, β, pi/2) && check_equal(γ, 2*pi/3)
end

function isRhombohedral(a::Real, b::Real, c::Real,
                       α::Real, β::Real, γ::Real)
    return check_equal(a, b, c) && check_equal(α, β, pi/2) &&\
           check_equal(γ, 2*pi/3)
end

function isOrthohombic(a::Real, b::Real, c::Real,
                       α::Real, β::Real, γ::Real)
    return check_not_equal(a, b, c) && check_equal(α, β, γ, pi/2)
end

function isMonoclinic(a::Real, b::Real, c::Real,
                      α::Real, β::Real, γ::Real)
    return check_not_equal(a, c) && check_equal(a, b) &&\
           check_equal(α, γ, pi/2) && check_not_equal(β, pi/2)
end

Base.Bool(c::Crystal) = true # For ease of testing

function volume(cl::Crystal)
    return (cl.a * cl.b * cl.c *
           sqrt( 1+2*cos(cl.α)*cos(cl.β)*cos(cl.γ) -
           cos(cl.α)^2 - cos(cl.β)^2 - cos(cl.γ)^2 ) )
end

function volume(cl::Monoclinic)
    cl.a * cl.b * cl.c * sin(cl.β)
end

function volume(cl::Union{Orthorhombic, Tetragonal, Cubic})
    cl.a * cl.b * cl.c
end

deg_to_rad(deg::Real) = deg/180*pi

function get_crystal(lattice_param::AbstractVector)
    length(lattice_param) == 6 || throw("There should be 6 lattice parameters!")
    a, b, c, α, β, γ = lattice_param
    α, β, γ = deg_to_rad.((α, β, γ))
    t = typeof(a)
    if isCubic(a, b, c, α, β, γ)
       return Cubic{t}(a)
    elseif isTetragonal(a, b, c, α, β, γ)
       return Tetragonal{t}(a, c)
    elseif isHexagonal(a, b, c, α, β, γ)
       return Hexagonal{t}(a, c)
    elseif isRhombohedral(a, b, c, α, β, γ)
       return Rhombohedral{t}(a) # TODO Or use trigonal????
    elseif isOrthohombic(a, b, c, α, β, γ)
       return Orthohombic{t}(a, b, c)
    elseif isMonoclinic(a, b, c, α, β, γ)
       return Monoclinic{t}(a, c, β)
    else
       return Triclinic{t}(a, b, c, α, β, γ)
    end
end

# (Crystal)(Peak) gives the peak position
# Fallback function and for triclinic
function (cl::Crystal)(P::Peak)
    (2pi/volume(cl) *
    sqrt(P.h^2 * cl.b^2 * cl.c^2 * sin(cl.α)^2
    + P.k^2 * cl.a^2 * cl.c^2 * sin(cl.β)^2
    + P.l^2 * cl.a^2 * cl.b^2 * sin(cl.γ)^2
    + 2*P.h * P.k * cl.a * cl.b * cl.c^2 * (cos(cl.α)*cos(cl.β) - cos(cl.γ))
    + 2*P.k * P.l * cl.a^2 * cl.b * cl.c * (cos(cl.β)*cos(cl.γ) - cos(cl.α))
    + 2*P.h * P.l * cl.a * cl.b^2 * cl.c * (cos(cl.α)*cos(cl.γ) - cos(cl.β))
    ))
end

function (cl::Cubic)(P::Peak)
    sqrt(P.h^2 + P.k^2 + P.l^2) / cl.a
end

# TODO More to go
