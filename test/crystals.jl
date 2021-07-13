using Test

include("../src/Crystal.jl")

@testset "Crystal object creation" begin
    @test_throws DomainError Triclinic{Float64}(1.0,1.0,3.0,4.0,5.0,6.0)
    @test_throws DomainError Triclinic{Float64}(1.0,2.0,3.0,4.0,5.0,5.0)
    @test Bool(Triclinic{Float64}(1.0,2.0,3.0,4.0,5.0,6.0))

    @test_throws DomainError Monoclinic{Float64}()
    @test_throws DomainError Monoclinic{Float64}()
    @test_throws DomainError Monoclinic{Float64}()
    @test_throws DomainError Monoclinic{Float64}()
    @test_throws DomainError Monoclinic{Float64}()
    @test Bool(Monoclinic{Float64}())

    @test_throws DomainError Orthorhombic{Float64}()
    @test_throws DomainError Orthorhombic{Float64}()
    @test Bool(Orthorhombic{Float64}())

    @test_throws DomainError Tetragonal{Float64}()
    @test_throws DomainError Tetragonal{Float64}()
    @test_throws DomainError Tetragonal{Float64}()
    @test_throws DomainError Tetragonal{Float64}()
    @test Bool(Tetragonal{Float64}())
end
