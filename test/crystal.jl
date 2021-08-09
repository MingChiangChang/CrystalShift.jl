using Test

include("../src/crystal.jl")

@testset "Utilities" begin
    @test check_not_equal(1,2,3)
    @test !check_not_equal(1,2,3,4,5,6,7,8,1)
    @test check_equal(1,1,1,1,1,1,1,1)
    @test !check_equal(1,1,1,1,1,1,2)

    # TODO isFuncs tests here
    # TODO deg_to_rad test
end

@testset "Crystal object creation" begin
    @test_throws MethodError Triclinic{Float64}(1.0,1.0,3.0,4.0,5.0,6.0)
    @test_throws MethodError Triclinic{Float64}(1.0,2.0,3.0,4.0,5.0,5.0)
    @test Bool(Triclinic{Float64}(1.0,2.0,3.0,4.0,5.0,6.0))

    @test_throws MethodError Monoclinic{Float64}(1.0, 1.0, 1.0)
    @test_throws MethodError Monoclinic{Float64}(2.0, 1.0, 1.0, pi/2, 1.0, pi/2)
    @test_throws MethodError Monoclinic{Float64}(1.0, 2.0, pi/2)
    @test Bool(Monoclinic{Float64}(1.0, 2.0, 1.0))

    @test_throws MethodError Orthorhombic{Float64}(1.0, 1.0, 2.0)
    @test_throws MethodError Orthorhombic{Float64}(1.0, 2.0, 1.0)
    @test_throws MethodError Orthorhombic{Float64}(1.0, 2.0, 1.0)
    @test_throws MethodError Orthorhombic{Float64}(2.0, 3.0, 1.0, 1/2, pi/2, pi/2)
    @test Bool(Orthorhombic{Float64}(2.0, 3.0, 1.0))

    @test_throws MethodError Tetragonal{Float64}(2.0, 2.0)
    @test_throws MethodError Tetragonal{Float64}(1.0, 2.0, 1.0, pi/2, 1/2, pi/2)
    @test Bool(Tetragonal{Float64}(1.0, 2.0))

    @test_throws MethodError Rhombohedral{Float64}(1.0, 1.0, 1.0, pi/2, pi/2, pi/3)
    @test Bool(Rhombohedral{Float64}(1.0))

    @test_throws MethodError Hexagonal{Float64}(2.0, 2.0, 2.0)
    @test_throws MethodError Hexagonal{Float64}(2.0, 2.0, 3.0, pi/2, pi/2, pi/3)
    @test Bool(Hexagonal{Float64}(2.0, 3.0))

    @test_throws MethodError Cubic{Float64}(1.0, 2.0, 1.0)
    @test_throws MethodError Cubic{Float64}(1.0, 1.0, 2.0)
    @test_throws MethodError Cubic{Float64}(2.0, 3.0, 1.0)
    @test_throws MethodError Cubic{Float64}(1.0, 1.0, 1.0, pi/2, pi/2, 1/2)
    @test Bool(Cubic{Float64}(1.0))

    # TODO get_crystal test
end

@testset "Volume" begin

end

@testset "peak position" begin

end
