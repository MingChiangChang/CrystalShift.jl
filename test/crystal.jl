using Test

include("../src/crystal.jl")

@testset "Utilities" begin
    @test check_not_equal(1,2,3)
    @test !check_not_equal(1,2,3,4,5,6,7,8,1)
    @test check_equal(1,1,1,1,1,1,1,1)
    @test !check_equal(1,1,1,1,1,1,2)
end

@testset "Crystal object creation" begin
    @test_throws DomainError Triclinic{Float64}(1.0,1.0,3.0,4.0,5.0,6.0)
    @test_throws DomainError Triclinic{Float64}(1.0,2.0,3.0,4.0,5.0,5.0)
    @test Bool(Triclinic{Float64}(1.0,2.0,3.0,4.0,5.0,6.0))

    @test_throws DomainError Monoclinic{Float64}(1.0, 1.0, 1.0, pi/2, 1.0, pi/2)
    @test_throws DomainError Monoclinic{Float64}(2.0, 1.0, 1.0, pi/2, 1.0, pi/2)
    @test_throws DomainError Monoclinic{Float64}(1.0, 3.0, 2.0, pi/2, 1.0, pi/2)
    @test_throws DomainError Monoclinic{Float64}(1.0, 1.0, 2.0, 1.0, 1.0, pi/2)
    @test_throws DomainError Monoclinic{Float64}(1.0, 1.0, 2.0, pi/2, 1.0, 0.1)
    @test_throws DomainError Monoclinic{Float64}(1.0, 1.0, 2.0, pi/2, pi/2, pi/2)
    @test Bool(Monoclinic{Float64}(1.0, 1.0, 2.0, pi/2, 1.0, pi/2))

    @test_throws DomainError Orthorhombic{Float64}(1.0, 1.0, 2.0, pi/2, pi/2, pi/2)
    @test_throws DomainError Orthorhombic{Float64}(1.0, 2.0, 1.0, pi/2, pi/2, pi/2)
    @test_throws DomainError Orthorhombic{Float64}(1.0, 2.0, 1.0, pi/2, pi/2, pi/2)
    @test_throws DomainError Orthorhombic{Float64}(2.0, 3.0, 1.0, 1/2, pi/2, pi/2)
    @test_throws DomainError Orthorhombic{Float64}(2.0, 3.0, 1.0, pi/2, 1/2, pi/2)
    @test_throws DomainError Orthorhombic{Float64}(2.0, 3.0, 1.0, pi/2, pi/2, 1/2)
    @test Bool(Orthorhombic{Float64}(2.0, 3.0, 1.0, pi/2, pi/2, pi/2))

    @test_throws DomainError Tetragonal{Float64}(1.0, 2.0, 1.0, pi/2, pi/2, pi/2)
    @test_throws DomainError Tetragonal{Float64}(1.0, 2.0, 2.0, pi/2, pi/2, pi/2)
    @test_throws DomainError Tetragonal{Float64}(1.0, 2.0, 3.0, pi/2, pi/2, pi/2)
    @test_throws DomainError Tetragonal{Float64}(1.0, 2.0, 1.0, pi/2, 1/2, pi/2)
    @test Bool(Tetragonal{Float64}(1.0, 1.0, 2.0, pi/2, pi/2, pi/2))

    @test_throws DomainError Rhombohedral{Float64}(1.0, 2.0, 1.0, pi/2, pi/2, 2*pi/3)
    @test_throws DomainError Rhombohedral{Float64}(1.0, 1.0, 2.0, pi/2, pi/2, 2*pi/3)
    @test_throws DomainError Rhombohedral{Float64}(2.0, 1.0, 1.0, pi/2, pi/2, 2*pi/3)
    @test_throws DomainError Rhombohedral{Float64}(1.0, 1.0, 1.0, 1/2, pi/2, 2*pi/3)
    @test_throws DomainError Rhombohedral{Float64}(1.0, 1.0, 1.0, pi/2, 1/2, 2*pi/3)
    @test_throws DomainError Rhombohedral{Float64}(1.0, 1.0, 1.0, pi/2, pi/2, pi/3)
    @test Bool(Rhombohedral{Float64}(1.0, 1.0, 1.0, pi/2, pi/2, 2*pi/3))

    @test_throws DomainError Hexagonal{Float64}(2.0, 2.0, 2.0, pi/2, pi/2, 2*pi/3)
    @test_throws DomainError Hexagonal{Float64}(1.0, 2.0, 3.0, pi/2, pi/2, 2*pi/3)
    @test_throws DomainError Hexagonal{Float64}(2.0, 2.0, 3.0, 1/2, pi/2, 2*pi/3)
    @test_throws DomainError Hexagonal{Float64}(2.0, 2.0, 3.0, pi/2, 1/2, 2*pi/3)
    @test_throws DomainError Hexagonal{Float64}(2.0, 2.0, 3.0, pi/2, pi/2, pi/3)
    @test Bool(Hexagonal{Float64}(2.0, 2.0, 3.0, pi/2, pi/2, 2*pi/3))

    @test_throws DomainError Cubic{Float64}(1.0, 2.0, 1.0, pi/2, pi/2, pi/2)
    @test_throws DomainError Cubic{Float64}(1.0, 1.0, 2.0, pi/2, pi/2, pi/2)
    @test_throws DomainError Cubic{Float64}(2.0, 3.0, 1.0, pi/2, pi/2, pi/2)
    @test_throws DomainError Cubic{Float64}(1.0, 1.0, 1.0, 1/2, pi/2, pi/2)
    @test_throws DomainError Cubic{Float64}(1.0, 1.0, 1.0, pi/2, 1/2, pi/2)
    @test_throws DomainError Cubic{Float64}(1.0, 1.0, 1.0, pi/2, pi/2, 1/2)
    @test Bool(Cubic{Float64}(1.0, 1.0, 1.0, pi/2, pi/2, pi/2))
end
