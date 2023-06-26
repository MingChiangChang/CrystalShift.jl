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

    sincos_α::Tuple{T, T}
    sincos_β::Tuple{T, T}
    sincos_γ::Tuple{T, T}

    volume::T

    free_param::Int8
    # Include sincos_alpha = sincos(α )
    # volume = volume(itself)
    function Triclinic{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:Real}
        #check_not_equal(a, b, c) ||  error("a,b,c should be different for triclinic")
        #check_not_equal(α, β, γ, pi/2) || error("α, β, γ, pi/2 should be different for triclinic")
        new{T}(a, b, c, α, β, γ, sincos(α), sincos(β), sincos(γ),
              volume(a, b, c, α, β, γ), 6)
    end
end

struct Monoclinic{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    sincos_α::Tuple{T, T}
    sincos_β::Tuple{T, T}
    sincos_γ::Tuple{T, T}

    volume::T

    free_param::Int8

    function Monoclinic{T}(a::T, b::T, c::T, β::T) where {T<:Real}
        #check_not_equal(a, b, c) || error("a,c should be different for monoclinic")
        #check_not_equal(β, pi/2) || error("β, pi/2 should be different for monoclinic")
        new{T}(a, b, c, pi/2, β, pi/2,
               (1., 0.), sincos(β), (1., 0.),
               a * b * c * sin(β), 4)
    end
end

struct Orthorhombic{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    sincos_α::Tuple{T, T}
    sincos_β::Tuple{T, T}
    sincos_γ::Tuple{T, T}

    volume::T

    free_param::Int8

    function Orthorhombic{T}(a::T, b::T, c::T) where {T<:Real}
        #check_not_equal(a, b, c) || error("a, b, c should be different for orthhombic")
        new{T}(a, b, c, pi/2, pi/2, pi/2, 
               (1., 0.), (1., 0.), (1., 0.), 
               a*b*c, 3)
    end
end

struct Tetragonal{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    sincos_α::Tuple{T, T}
    sincos_β::Tuple{T, T}
    sincos_γ::Tuple{T, T}

    volume::T

    free_param::Int8

    function Tetragonal{T}(a::T, c::T) where {T<:Real}
        #check_not_equal(a, c) || error("a, c should be different for tetragonal")
        new{T}(a, a, c, pi/2, pi/2, pi/2, 
               (1., 0.), (1., 0.), (1., 0.),
               a*a*c, 2)
    end
end

struct Rhombohedral{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    sincos_α::Tuple{T, T}
    sincos_β::Tuple{T, T}
    sincos_γ::Tuple{T, T}

    volume::T

    free_param::Int8

    function Rhombohedral{T}(a, α) where {T<:Real}
        sincos_α = sincos(α)
        new{T}(a, a, a, α, α, α, 
        sincos_α, sincos_α, sincos_α,
               volume(a, a, a, α, α, α), 2)
    end
end

struct Hexagonal{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    sincos_α::Tuple{T, T}
    sincos_β::Tuple{T, T}
    sincos_γ::Tuple{T, T}

    volume::T

    free_param::Int8

    function Hexagonal{T}(a::T, c::T) where {T<:Real}
        #check_not_equal(a, c) || error("a, c should be different for hexagonal")
        new{T}(a, a, c, pi/2, pi/2, 2*pi/3,
               (1., 0.), (1., 0.), sincos(2*pi/3),
               volume(a, a, c, pi/2, pi/2, 2*pi/3), 2)
    end
end

struct Cubic{T}<:Crystal{T}
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    sincos_α::Tuple{T, T}
    sincos_β::Tuple{T, T}
    sincos_γ::Tuple{T, T}

    volume::T

    free_param::Int8

    # get_property to get b,c, α, β, γ
    function Cubic{T}(a::T) where {T<:Real}
        new{T}(a, a, a, pi/2, pi/2, pi/2, 
               (1., 0.), (1., 0.), (1., 0.),
               a*a*a, 1)
    end
end

function Base.show(io::IO, ct::Crystal)
    println("a: $(ct.a), b: $(ct.b), c: $(ct.c)")
    println("α: $(ct.α/pi*180), β: $(ct.β/pi*180), γ: $(ct.γ/pi*180)")
end

check_not_equal(x...) = length(Set(x)) == length(x)
check_equal(x...) = all(y->isapprox(y, x[1], rtol=0.001), x)

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

volume(cl::Crystal) = volume(cl.a, cl.b, cl.c, cl.α, cl.β, cl.γ)

function volume(cl::Monoclinic)
    cl.a * cl.b * cl.c * cl.sincos_β[1]
end

function volume(cl::Union{Orthorhombic, Tetragonal, Cubic})
    cl.a * cl.b * cl.c
end

# deg_to_rad(deg::Real) = deg/180*pi

function get_crystal(lattice_param::AbstractVector, deg::Bool=true)
    length(lattice_param) == 6 || error("There should be 6 lattice parameters!")
    a, b, c, α, β, γ = lattice_param

    if deg
        α, β, γ = deg2rad.((α, β, γ))
    end

    t = typeof(a)
    if isCubic(a, b, c, α, β, γ)
       return Cubic{t}(a)
    elseif isTetragonal(a, b, c, α, β, γ)
       return Tetragonal{t}(a, c)
    elseif isHexagonal(a, b, c, α, β, γ)
       return Hexagonal{t}(a, c)
    elseif isRhombohedral(a, b, c, α, β, γ)
       return Rhombohedral{t}(a, α)
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
    try
        (2pi/cl.volume *
        sqrt(P.h^2 * cl.b^2 * cl.c^2 * cl.sincos_α[1]^2
        + P.k^2 * cl.a^2 * cl.c^2 * cl.sincos_β[1]^2
        + P.l^2 * cl.a^2 * cl.b^2 * cl.sincos_γ[1]^2
        + 2*P.h * P.k * cl.a * cl.b * cl.c^2 * (cl.sincos_α[2]*cl.sincos_β[2] - cl.sincos_γ[2])
        + 2*P.k * P.l * cl.a^2 * cl.b * cl.c * (cl.sincos_β[2]*cl.sincos_γ[2] - cl.sincos_α[2])
        + 2*P.h * P.l * cl.a * cl.b^2 * cl.c * (cl.sincos_α[2]*cl.sincos_γ[2] - cl.sincos_β[2])
        ))
    catch DomainError
        return Inf
    end
end

function (cl::Cubic)(P::Peak)
    2pi*sqrt(P.h^2 + P.k^2 + P.l^2) / cl.a
end

function (cl::Tetragonal)(P::Peak)
    2pi*sqrt(P.h^2*cl.c^2 + P.k^2*cl.c^2 + P.l^2*cl.a^2) / (cl.a*cl.c)
end
function (cl::Hexagonal)(P::Peak)
    (2pi/cl.volume *
    sqrt(P.h^2 * cl.b^2 * cl.c^2
    + P.k^2 * cl.a^2 * cl.c^2
    + P.l^2 * cl.a^2 * cl.b^2 * cl.sincos_γ[1]^2
    + 2*P.h * P.k * cl.a * cl.b * cl.c^2 * (cl.sincos_α[2]*cl.sincos_β[2] - cl.sincos_γ[2])))
end

function (cl::Orthorhombic)(P::Peak)
    return (2pi/(cl.a*cl.b*cl.c) *
            sqrt(P.h^2 * cl.b^2 * cl.c^2
                + P.k^2 * cl.a^2 * cl.c^2
                + P.l^2 * cl.a^2 * cl.b^2))
end

function (cl::Rhombohedral)(P::Peak)
    (2pi/cl.volume *
    sqrt(P.h^2 * cl.a^4 * cl.sincos_α[1]^2
    + P.k^2 * cl.a^4 * cl.sincos_α[1]^2
    + P.l^2 * cl.a^4 * cl.sincos_α[1]^2
    + ((2cl.a^4)*(cl.sincos_α[2]^2 - cl.sincos_α[2])
    * (P.h * P.k + P.k * P.l + P.h * P.l) )))
end

function (cl::Monoclinic)(P::Peak)
    (2pi/cl.volume *
    sqrt(P.h^2 * cl.b^2 * cl.c^2
    + P.k^2 * cl.a^2 * cl.c^2 * cl.sincos_β[1]^2
    + P.l^2 * cl.a^2 * cl.b^2
    + 2*P.h * P.l * cl.a * cl.b^2 * cl.c * (- cl.sincos_β[2])))
end

get_free_lattice_params(cl::Cubic) = [cl.a]
get_free_lattice_params(cl::Tetragonal) = [cl.a, cl.c]
get_free_lattice_params(cl::Hexagonal) = [cl.a, cl.c]
get_free_lattice_params(cl::Orthorhombic) = [cl.a, cl.b, cl.c]
get_free_lattice_params(cl::Rhombohedral) = [cl.a, cl.α]
get_free_lattice_params(cl::Monoclinic) = [cl.a, cl.b, cl.c, cl.β]
get_free_lattice_params(cl::Triclinic) = [cl.a, cl.b, cl.c, cl.α, cl.β, cl.γ]

## Abstract type to serve as supertype for 7 different crystal systems
#abstract type Crystal{T} end
#
## Unit in angstrom and radians
## TODO use Base.getproperty to make things easier
## TODO: Include precalculation factors to improve speed
#struct Triclinic{T}<:Crystal{T}
#    a::T
#    b::T
#    c::T
#
#    α::T
#    β::T
#    γ::T
#
#    sincos_α::Tuple{T, T}
#    sincos_β::Tuple{T, T}
#    sincos_γ::Tuple{T, T}
#
#    # Precalculate
#    b_square_c_square_sin_α_square::T
#    a_square_c_square_sin_β_square::T
#    a_square_b_square_sin_γ_square::T
#    a_b_c_square_cos_α_cos_β_minus_cos_γ::T
#    b_c_a_square_cos_β_cos_γ_minus_cos_α::T
#    a_c_b_square_cos_α_cos_γ_minus_cos_β::T
#    volume::T
#
#    free_param::Int8
#    # Include sincos_alpha = sincos(α )
#    # volume = volume(itself)
#    function Triclinic{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:Real}
#        #check_not_equal(a, b, c) ||  error("a,b,c should be different for triclinic")
#        #check_not_equal(α, β, γ, pi/2) || error("α, β, γ, pi/2 should be different for triclinic")
#        new{T}(a, b, c, α, β, γ, sincos(α), sincos(β), sincos(γ),
#              (b*c*sin(α))^2, (a*c*sin(β))^2, (a*b*sin(γ))^2,
#              a*b*c^2*(cos(α)*cos(β)-cos(γ)),
#              b*c*a^2*(cos(β)*cos(γ)-cos(α)),
#              a*c*b^2*(cos(α)*cos(γ)-cos(β)),
#              volume(a, b, c, α, β, γ), 6)
#    end
#end
#
#struct Monoclinic{T}<:Crystal{T}
#    a::T
#    b::T
#    c::T
#
#    α::T
#    β::T
#    γ::T
#
#    # Precalculate
#    b_square_c_square::T
#    a_square_c_square_sin_β_square::T
#    a_square_b_square::T
#    a_b_square_c_cos_β::T
#
#    volume::T
#
#    free_param::Int8
#
#    function Monoclinic{T}(a::T, b::T, c::T, β::T) where {T<:Real}
#        #check_not_equal(a, b, c) || error("a,c should be different for monoclinic")
#        #check_not_equal(β, pi/2) || error("β, pi/2 should be different for monoclinic")
#        new{T}(a, b, c, pi/2, β, pi/2,
#               (b*c)^2, (a*c*sin(β))^2, (a*b)^2, a*b^2*c*cos(β),
#               a * b * c * sin(β), 4)
#    end
#end
#
#struct Orthorhombic{T}<:Crystal{T}
#    a::T
#    b::T
#    c::T
#
#    α::T
#    β::T
#    γ::T
#
#    # Precalculate
#    b_square_c_square::T
#    a_square_c_square::T
#    a_square_b_square::T
#    volume::T
#
#    free_param::Int8
#
#    function Orthorhombic{T}(a::T, b::T, c::T) where {T<:Real}
#        #check_not_equal(a, b, c) || error("a, b, c should be different for orthhombic")
#        new{T}(a, b, c, pi/2, pi/2, pi/2,
#               b^2*c^2, a^2*c^2, a^2*b^2,
#               a*b*c, 3)
#    end
#end
#
#struct Tetragonal{T}<:Crystal{T}
#    a::T
#    b::T
#    c::T
#
#    α::T
#    β::T
#    γ::T
#
#    # Precalculate
#    a_square::T
#    c_square::T
#    ac::T
#    volume::T
#
#    free_param::Int8
#
#    function Tetragonal{T}(a::T, c::T) where {T<:Real}
#        #check_not_equal(a, c) || error("a, c should be different for tetragonal")
#        new{T}(a, a, c, pi/2, pi/2, pi/2,
#               a^2, c^2, a*c, a*a*c, 2)
#    end
#end
#
#struct Rhombohedral{T}<:Crystal{T}
#    a::T
#    b::T
#    c::T
#
#    α::T
#    β::T
#    γ::T
#
#    # Precalculate
#    a_quad_sin_α_square::T
#    two_a_quad_cos_α_square_minus_cos_α::T
#    volume::T
#
#    free_param::Int8
#
#    function Rhombohedral{T}(a, α) where {T<:Real}
#        sincos_α = sincos(α)
#        new{T}(a, a, a, α, α, α,
#        a^4*sincos_α[1]^2,
#        2a^4*(sincos_α[2]^2 - sincos_α[2]),
#        volume(a, a, a, α, α, α), 2)
#    end
#end
#
#struct Hexagonal{T}<:Crystal{T}
#    a::T
#    b::T
#    c::T
#
#    α::T
#    β::T
#    γ::T
#
#    # Precalculate terms
#    a_square_c_square::T
#    a_quad_sin_γ_square::T
#    #ac::T
#    volume::T
#
#    free_param::Int8
#
#    function Hexagonal{T}(a::T, c::T) where {T<:Real}
#        #check_not_equal(a, c) || error("a, c should be different for hexagonal")
#        new{T}(a, a, c, pi/2, pi/2, 2*pi/3,
#               a^2*c^2, a^4*0.75,
#               volume(a, a, c, pi/2, pi/2, 2*pi/3), 2)
#    end
#end
#
#struct Cubic{T}<:Crystal{T}
#    a::T
#    b::T
#    c::T
#
#    α::T
#    β::T
#    γ::T
#
#    volume::T
#
#    free_param::Int8
#
#    # get_property to get b,c, α, β, γ
#    function Cubic{T}(a::T) where {T<:Real}
#        new{T}(a, a, a, pi/2, pi/2, pi/2,
#               a^3, 1)
#    end
#end
#
#function Base.show(io::IO, ct::Crystal)
#    println("a: $(ct.a), b: $(ct.b), c: $(ct.c)")
#    println("α: $(ct.α/pi*180), β: $(ct.β/pi*180), γ: $(ct.γ/pi*180)")
#end
#
#check_not_equal(x...) = length(Set(x)) == length(x)
#check_equal(x...) = all(y->isapprox(y, x[1], rtol=0.001), x)
#
#function isCubic(a::Real, b::Real, c::Real,
#                 α::Real, β::Real, γ::Real)
#    return check_equal(a, b, c) && check_equal(α, β, γ, pi/2)
#end
#
#function isTetragonal(a::Real, b::Real, c::Real,
#                      α::Real, β::Real, γ::Real)
#    return (check_not_equal(a, c) && check_equal(a, b) &&
#           check_equal(a, b) && check_equal(α, β, γ, pi/2))
#end
#
#function isHexagonal(a::Real, b::Real, c::Real,
#                      α::Real, β::Real, γ::Real)
#    return (check_not_equal(a, c) && check_equal(a, b) &&
#           check_equal(α, β, pi/2) && check_equal(γ, 2*pi/3))
#end
#
#function isRhombohedral(a::Real, b::Real, c::Real,
#                       α::Real, β::Real, γ::Real)
#    return check_equal(a, b, c) && check_equal(α, β, γ) &&
#           check_not_equal(α, pi/2)
#end
#
#function isOrthohombic(a::Real, b::Real, c::Real,
#                       α::Real, β::Real, γ::Real)
#    return check_not_equal(a, b, c) && check_equal(α, β, γ, pi/2)
#end
#
#function isMonoclinic(a::Real, b::Real, c::Real,
#                      α::Real, β::Real, γ::Real)
#    return (check_not_equal(a, b, c) &&
#           check_equal(α, γ, pi/2) && check_not_equal(β, pi/2))
#end
#
#Base.Bool(c::Crystal) = true # For ease of testing
#
#volume(cl::Crystal) = volume(cl.a, cl.b, cl.c, cl.α, cl.β, cl.γ)
#
#function volume(cl::Monoclinic)
#    cl.a * cl.b * cl.c * sin(cl.β)
#end
#
#function volume(cl::Union{Orthorhombic, Tetragonal, Cubic})
#    cl.a * cl.b * cl.c
#end
#
## deg_to_rad(deg::Real) = deg/180*pi
#
#function get_crystal(lattice_param::AbstractVector, deg::Bool=true)
#    length(lattice_param) == 6 || error("There should be 6 lattice parameters!")
#    a, b, c, α, β, γ = lattice_param
#
#    if deg
#        α, β, γ = deg2rad.((α, β, γ))
#    end
#
#    t = typeof(a)
#    if isCubic(a, b, c, α, β, γ)
#       return Cubic{t}(a)
#    elseif isTetragonal(a, b, c, α, β, γ)
#       return Tetragonal{t}(a, c)
#    elseif isHexagonal(a, b, c, α, β, γ)
#       return Hexagonal{t}(a, c)
#    elseif isRhombohedral(a, b, c, α, β, γ)
#       return Rhombohedral{t}(a, α)
#    elseif isOrthohombic(a, b, c, α, β, γ)
#       return Orthorhombic{t}(a, b, c)
#    elseif isMonoclinic(a, b, c, α, β, γ)
#       return Monoclinic{t}(a, b, c, β)
#    else
#       return Triclinic{t}(a, b, c, α, β, γ)
#    end
#end
#
## (Crystal)(Peak) gives the peak position
## Fallback function and for triclinic
#function (cl::Crystal)(P::Peak)
#    try
#        (2pi/cl.volume *
#        sqrt(P.h^2 * cl.b_square_c_square_sin_α_square
#        + P.k^2 * cl.a_square_c_square_sin_β_square
#        + P.l^2 * cl.a_square_b_square_sin_γ_square
#        + 2*P.h * P.k * cl.a_b_c_square_cos_α_cos_β_minus_cos_γ
#        + 2*P.k * P.l * cl.b_c_a_square_cos_β_cos_γ_minus_cos_α
#        + 2*P.h * P.l * cl.a_c_b_square_cos_α_cos_γ_minus_cos_β
#        ))
#    catch DomainError
#        return Inf
#    end
#end
#
#function (cl::Cubic)(P::Peak)
#    2pi*sqrt(P.h^2 + P.k^2 + P.l^2) / cl.a
#end
#
#function (cl::Tetragonal)(P::Peak)
#    2pi*sqrt(P.h^2*cl.c_square + P.k^2*cl.c_square + P.l^2*cl.a_square) / (cl.ac)
#end
#function (cl::Hexagonal)(P::Peak)
#    (2pi/cl.volume *
#    sqrt(P.h^2 * cl.a_square_c_square
#    + P.k^2 * cl.a_square_c_square
#    + P.l^2 * cl.a_quad_sin_γ_square
#    + P.h * P.k * cl.a_square_c_square))
#end
#
#function (cl::Orthorhombic)(P::Peak)
#    # try
#        return (2pi/cl.volume *
#                sqrt(P.h^2 * cl.b_square_c_square
#                    + P.k^2 * cl.a_square_c_square
#                    + P.l^2 * cl.a_square_b_square))
#    # catch DomainError
#    #     return Inf
#    # end
#end
#
#function (cl::Rhombohedral)(P::Peak)
#    (2pi/cl.volume *
#    sqrt( (P.h^2 + P.k^2 + P.l^2) * cl.a_quad_sin_α_square
#    + (cl.two_a_quad_cos_α_square_minus_cos_α
#    * (P.h * P.k + P.k * P.l + P.h * P.l) )))
#end
#
#@fastmath function (cl::Monoclinic)(P::Peak)
#    (2pi/cl.volume *
#    sqrt(P.h^2 * cl.b_square_c_square
#    + P.k^2 * cl.a_square_c_square_sin_β_square
#    + P.l^2 * cl.a_square_b_square
#    - 2*P.h * P.l * cl.a_b_square_c_cos_β))
#    #+ 2*P.h * P.l * cl.a * cl.b^2 * cl.c * (- cl.sincos_β[2])))
#end
#
#get_free_lattice_params(cl::Cubic) = [cl.a]
#get_free_lattice_params(cl::Tetragonal) = [cl.a, cl.c]
#get_free_lattice_params(cl::Hexagonal) = [cl.a, cl.c]
#get_free_lattice_params(cl::Orthorhombic) = [cl.a, cl.b, cl.c]
#get_free_lattice_params(cl::Rhombohedral) = [cl.a, cl.α]
#get_free_lattice_params(cl::Monoclinic) = [cl.a, cl.b, cl.c, cl.β]
#get_free_lattice_params(cl::Triclinic) = [cl.a, cl.b, cl.c, cl.α, cl.β, cl.γ]
