# Abstract type to serve as supertype for 7 different crystal systems
abstract type Crystal{T} end

# Unit in angstrom and radians
# TODO use Base.getproperty to make things easier
struct Triclinic{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param::Int8

    function Triclinic{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:Real}
        check_not_equal(a, b, c) ||  error("a,b,c should be different for triclinic")
        check_not_equal(α, β, γ, pi/2) || error("α, β, γ, pi/2 should be different for triclinic")
        new{T}(a, b, c, α, β, γ, 6)
    end
end

struct Monoclinic{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param::Int8

    function Monoclinic{T}(a::T, b::T, c::T, β::T) where {T<:Real}
        check_not_equal(a, b, c) || error("a,c should be different for monoclinic")
        check_not_equal(β, pi/2) || error("β, pi/2 should be different for monoclinic")
        new{T}(a, b, c, pi/2, β, pi/2, 4)
    end
end

struct Orthorhombic{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param::Int8

    function Orthorhombic{T}(a::T, b::T, c::T) where {T<:Real}
        check_not_equal(a, b, c) || error("a, b, c should be different for orthhombic")
        new{T}(a, b, c, pi/2, pi/2, pi/2, 3)
    end
end

struct Tetragonal{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param::Int8

    function Tetragonal{T}(a::T, c::T) where {T<:Real}
        check_not_equal(a, c) || error("a, c should be different for tetragonal")
        new{T}(a, a, c, pi/2, pi/2, pi/2, 2)
    end
end

struct Rhombohedral{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param::Int8

    function Rhombohedral{T}(a, α) where {T<:Real}
        new{T}(a, a, a, α, α, α, 2)
    end
end

struct Hexagonal{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param::Int8

    function Hexagonal{T}(a::T, c::T) where {T<:Real}
        check_not_equal(a, c) || error("a, c should be different for hexagonal")
        new{T}(a, a, c, pi/2, pi/2, 2*pi/3, 2)
    end
end

struct Cubic{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    free_param::Int8

    # get_property to get b,c, α, β, γ
    function Cubic{T}(a::T) where {T<:Real}
        new{T}(a, a, a, pi/2, pi/2, pi/2, 1)
    end
end

function Base.show(io::IO, ct::Crystal)
    println("a: $(ct.a), b: $(ct.b), c: $(ct.c)")
    println("α: $(ct.α/pi*180), β: $(ct.β/pi*180), γ: $(ct.γ/pi*180)")
end

check_not_equal(x...) = length(Set(x)) == length(x)
check_equal(x...) = all(y->y==x[1], x)

function isCubic(a::Real, b::Real, c::Real,
                 α::Real, β::Real, γ::Real)
    return check_equal(a, b, c) && check_equal(α, β, γ, pi/2)
end

function isTetragonal(a::Real, b::Real, c::Real,
                      α::Real, β::Real, γ::Real)
    return (check_not_equal(a, c) && check_equal(a, b) &&
           check_equal(a, b) && check_equal(α, β, γ, pi/2))
end

function isHexagonal(a::Real, b::Real, c::Real,
                      α::Real, β::Real, γ::Real)
    return (check_not_equal(a, c) && check_equal(a, b) &&
           check_equal(α, β, pi/2) && check_equal(γ, 2*pi/3))
end

function isRhombohedral(a::Real, b::Real, c::Real,
                       α::Real, β::Real, γ::Real)
    return check_equal(a, b, c) && check_equal(α, β, γ) &&
           check_not_equal(α, pi/2)
end

function isOrthohombic(a::Real, b::Real, c::Real,
                       α::Real, β::Real, γ::Real)
    return check_not_equal(a, b, c) && check_equal(α, β, γ, pi/2)
end

function isMonoclinic(a::Real, b::Real, c::Real,
                      α::Real, β::Real, γ::Real)
    return (check_not_equal(a, b, c) &&
           check_equal(α, γ, pi/2) && check_not_equal(β, pi/2))
end

Base.Bool(c::Crystal) = true # For ease of testing

function volume(cl::Crystal)
    return (cl.a * cl.b * cl.c
           * sqrt( 1+2*cos(cl.α)*cos(cl.β)*cos(cl.γ)
           - cos(cl.α)^2 - cos(cl.β)^2 - cos(cl.γ)^2 ) )
end

function volume(cl::Monoclinic)
    cl.a * cl.b * cl.c * sin(cl.β)
end

function volume(cl::Union{Orthorhombic, Tetragonal, Cubic})
    cl.a * cl.b * cl.c
end

deg_to_rad(deg::Real) = deg/180*pi

function get_crystal(lattice_param::AbstractVector, deg::Bool=true)
    length(lattice_param) == 6 || error("There should be 6 lattice parameters!")
    a, b, c, α, β, γ = lattice_param
    if deg
        α, β, γ = deg_to_rad.((α, β, γ))
    end
    t = typeof(a)
    if isCubic(a, b, c, α, β, γ)
       return Cubic{t}(a)
    elseif isTetragonal(a, b, c, α, β, γ)
       return Tetragonal{t}(a, c)
    elseif isHexagonal(a, b, c, α, β, γ)
       return Hexagonal{t}(a, c)
    elseif isRhombohedral(a, b, c, α, β, γ)
       return Rhombohedral{t}(a, c)
    elseif isOrthohombic(a, b, c, α, β, γ)
       return Orthorhombic{t}(a, b, c)
    elseif isMonoclinic(a, b, c, α, β, γ)
       return Monoclinic{t}(a, b, c, β)
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
    2pi*sqrt(P.h^2 + P.k^2 + P.l^2) / cl.a
end

function (cl::Tetragonal)(P::Peak)
    2pi*sqrt(P.h^2*cl.c^2 + P.k^2*cl.c^2 + P.l^2*cl.a^2) / (cl.a*cl.c)
end
function (cl::Hexagonal)(P::Peak)
    (2pi/volume(cl) *
    sqrt(P.h^2 * cl.b^2 * cl.c^2
    + P.k^2 * cl.a^2 * cl.c^2
    + P.l^2 * cl.a^2 * cl.b^2 * sin(cl.γ)^2
    + 2*P.h * P.k * cl.a * cl.b * cl.c^2 * (cos(cl.α)*cos(cl.β) - cos(cl.γ))))
end

function (cl::Orthorhombic)(P::Peak)
    (2pi/(cl.a*cl.b*cl.c) *
    sqrt(P.h^2 * cl.b^2 * cl.c^2
         + P.k^2 * cl.a^2 * cl.c^2
         + P.l^2 * cl.a^2 * cl.b^2))
end

function (cl::Rhombohedral)(P::Peak)
    (2pi/volume(cl) *
    sqrt(P.h^2 * cl.a^4 * sin(cl.α)^2
    + P.k^2 * cl.a^4 * sin(cl.α)^2
    + P.l^2 * cl.a^4 * sin(cl.α)^2
    + ((2cl.a^4)(cos(cl.α)^2 - cos(cl.α))
    * (P.h * P.k + P.k * P.l + P.h * P.l) )))
end

function (cl::Monoclinic)(P::Peak)
    (2pi/volume(cl) *
    sqrt(P.h^2 * cl.b^2 * cl.c^2
    + P.k^2 * cl.a^2 * cl.c^2 * sin(cl.β)^2
    + P.l^2 * cl.a^2 * cl.b^2))
end

get_free_params(cl::Cubic) = [cl.a]
get_free_params(cl::Tetragonal) = [cl.a, cl.c]
get_free_params(cl::Hexagonal) = [cl.a, cl.c]
get_free_params(cl::Orthorhombic) = [cl.a, cl.b, cl.c]
get_free_params(cl::Rhombohedral) = [cl.a, cl.α]
get_free_params(cl::Monoclinic) = [cl.a, cl.b, cl.c, cl.β]
get_free_params(cl::Triclinic) = [cl.a, cl.b, cl.c, cl.α, cl.β, cl.γ]