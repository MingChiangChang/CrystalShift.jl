using Test

include("../src/Crystal.jl")

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
end
