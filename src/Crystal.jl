using PhaseMapping

# Abstract type to as super type for 7 different
abstract type Crystal end

# TODO Constructors must make sure things that should be different
#      are actually different
struct Triclinic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    function Triclinic(a, b, c, α, β, γ)
        check_not_equal()
        new{}
    end
end

struct Monoclinic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
end

struct Orthorhombic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
end

struct Tetragonal{T}<:Crystal
    a::T
    c::T

    α::T
end

struct Trigonal{T}<:Crystal
    a::T
    c::T

    α::T
    β::T
end

struct Hexagonal{T}<:Crystal
    a::T
    c::T

    α::T
    β::T
end

struct Cubic{T}<:Crystal
    a::T

    α::T
end

function check_not_equal()

end
