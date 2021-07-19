using PhaseMapping

# Abstract type to as supertype for 7 different crystal systems
abstract type Crystal end

struct Triclinic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    function Triclinic{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:AbstractFloat}
        check_not_equal(a, b, c) || throw(DomainError((a, b, c), "a, b and c must not equal for Triclinic crystals")) #error("a, b and c must not equal for Triclinic crystals.")
        check_not_equal(α, β, γ) || throw(DomainError((α, β, γ), "α, β and γ must not equal for Triclinic crystals."))
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

    function Monoclinic{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:AbstractFloat}
        check_not_equal(a, c) || throw(DomainError((a, c), "a and c must not equal for monoclinic crystals."))
        check_not_equal(b, c) || throw(DomainError((b, c), "b and c must not equal for monoclinic crystals."))
        check_equal(a, b) || throw(DomainError((a, b), "a and b must be equal for monoclinic crystals"))
        check_equal(α, γ, pi/2) || throw(DomainError((α, γ, pi/2), "α, β and γ must not equal for monoclinic crystals."))
        check_not_equal(β, pi/2) || throw(DomainError((β, pi/2), "α, β and γ must not equal for monoclinic crystals."))
        new{T}(a, b, c, α, β, γ)
    end
end

struct Orthorhombic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    function Orthorhombic{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:AbstractFloat}
        check_not_equal(a, b, c) || throw(DomainError((a, b, c), "a and c must not equal for orthorhombic crystals."))
        check_equal(α, β, γ, pi/2) || throw(DomainError((α, β, γ, pi/2), "α, β and γ must equal π/2 for orthorhombic crystals."))
        new{T}(a, b, c, α, β, γ)
    end
end

struct Tetragonal{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    function Tetragonal{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:AbstractFloat}
        check_not_equal(a, c) || throw(DomainError((a, c), "a and c must not equal for tetragonal crystals."))
        check_not_equal(b, c) || throw(DomainError((b, c), "b and c must not equal for tetragonal crystals."))
        check_equal(a, b) || throw(DomainError((a, b), "a and b must be equal for tetragonal crystals."))
        check_equal(α, β, γ, pi/2) || throw(DomainError((α, β, γ, pi/2), "α, β and γ must equal π/2 for tetragonal crystals."))
        new{T}(a, b, c, α, β, γ)
    end
end

struct Rhombohedral{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    function Rhombohedral{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:AbstractFloat}
        check_equal(a, b, c) || throw(DomainError((a, b, c), "a, b and c must equal for rhombohedral crystals."))
        check_equal(α, β, pi/2) || throw(DomainError((α, β, pi/2), "α and β must equal π/2 for rhombohedral crystals."))
        check_equal(γ, 2*pi/3) || throw(DomainError((γ, 2*pi/3), "γ must equal to 2/3π for rhombohedral crystals"))
        new{T}(a, b, c, α, β, γ)
    end
end

struct Hexagonal{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    function Hexagonal{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:AbstractFloat}
        check_not_equal(a, c) || throw(DomainError((a, c), "a and c must not equal for hexagonal crystals."))
        check_not_equal(b, c) || throw(DomainError((b, c), "b and c must not equal for hexagonal crystals."))
        check_equal(a, b) || throw(DomainError((a,b), "a and b must be equal for hexagonal crystals."))
        check_equal(α, β, pi/2) || throw(DomainError((α, β, pi/2), "α and β must equal π/2 for hexagonal crystals."))
        check_equal(γ, 2*pi/3) || throw(DomainError((γ, 2*pi/3), "γ must eqaul 2/3π for hexagonal crystals."))
        new{T}(a, b, c, α, β, γ)
    end
end

struct Cubic{T}<:Crystal
    a::T
    b::T
    c::T

    α::T
    β::T
    γ::T

    function Cubic{T}(a::T, b::T, c::T, α::T, β::T, γ::T) where {T<:AbstractFloat}
        check_equal(a, b, c) || throw(DomainError((a, b, c), "a, b and c must equal for cubic crystals."))
        check_equal(α, β, γ, pi/2) || throw(DomainError((α, β, γ, pi/2), "α, β and γ must equal π/2 for cubic crystals"))
        new{T}(a, b, c, α, β, γ)
    end
end

check_not_equal(x...) = length(Set(x)) == length(x)
check_equal(x...) = all(y->y==x[1], x)

Cubic(a::AbstractFloat) = Cubic(a, a, a, pi/2, pi/2, pi/2)

Base.Bool(c::Crystal) = true # For ease of testing

function volume(cl::Crystal)
    cl.a*cl.b*cl.c*sqrt(1+2*cos(cl.α)*cos(cl.β)*cos(cl.γ)-cos(cl.α)^2-cos(cl.β)^2-cos(cl.γ)^2)
end

function volume(cl::Monoclinic)
    cl.a*cl.b*cl.c*sin(cl.β)
end

function volume(cl::Union{Orthorhombic, Tetragonal, Cubic})
    cl.a*cl.b*cl.c
end
