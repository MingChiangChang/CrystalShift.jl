using Test
using CrystalShift
using CrystalShift: check_equal, check_not_equal
using CrystalShift: get_crystal, volume, twoθ2q, get_free_params

@testset "Utilities" begin
    @test check_not_equal(1,2,3)
    @test !check_not_equal(1,2,3,4,5,6,7,8,1)
    @test check_equal(1,1,1,1,1,1,1,1)
    @test !check_equal(1,1,1,1,1,1,2)

    # TODO isFuncs tests here
    # @test deg_to_rad(90.) == pi/2

    @test isCubic(1, 1, 1, pi/2, pi/2, pi/2)
    @test !isCubic(1, 2, 1, pi/2, pi/2, pi/2)

    @test isTetragonal(1, 1, 2, pi/2, pi/2, pi/2)
    @test !isTetragonal(1, 2, 2, pi/2, pi/2, pi/2)

    @test isHexagonal(1,1,3, pi/2, pi/2, 2pi/3)
    @test !isHexagonal(1,2,3, pi/2, pi/2, 2pi/3)

    @test isRhombohedral(1, 1, 1, 0.5, 0.5, 0.5)
    @test !isRhombohedral(1, 1, 1, pi/2, pi/2, pi/2)
    
    @test isOrthohombic(1, 2, 3, pi/2, pi/2, pi/2)
    @test !isOrthohombic(1, 1, 1, pi/2, pi/2, pi/2)
    @test !isOrthohombic(1, 2, 3, 0.5, pi/2, pi/2)

    @test isMonoclinic(1, 2, 3, pi/2, 0.5, pi/2)
    @test !isMonoclinic(1, 2, 3, 0.5, pi/2, pi/2)
    @test !isMonoclinic(1, 2, 2, pi/2, 0.5, pi/2)
end

@testset "Volume and get_crystal" begin
    cl = get_crystal([12.45245, 3.08297, 5.87615, 90.000, 103.6836, 90.0000], true)
    @test typeof(cl) <: Monoclinic{Float64}
    @test isapprox(volume(cl), 219.185519, rtol=0.01)

    cl = get_crystal([5.39991, 5.39991, 5.39991, 55.8716, 55.8716, 55.8716], true)
    @test typeof(cl) <: Rhombohedral{Float64}
    @test isapprox(volume(cl), 100.683545, rtol=0.01)

    cl = get_crystal([5., 5., 5., 90., 90., 90.], true)
    @test typeof(cl) <: Cubic{Float64}
    @test volume(cl) == 125.
end

@testset "Crystal object creation" begin
    @test Bool(Triclinic{Float64}(1.0,2.0,3.0,4.0,5.0,6.0))

    @test_throws MethodError Monoclinic{Float64}(2.0, 1.0, 1.0, pi/2, 1.0, pi/2)
    @test Bool(Monoclinic{Float64}(1.0, 1.5, 2.0, 1.0))

    @test_throws MethodError Orthorhombic{Float64}(2.0, 3.0, 1.0, 1/2, pi/2, pi/2)
    @test Bool(Orthorhombic{Float64}(2.0, 3.0, 1.0))

    @test_throws MethodError Tetragonal{Float64}(1.0, 2.0, 1.0, pi/2, 1/2, pi/2)
    @test Bool(Tetragonal{Float64}(1.0, 2.0))

    @test_throws MethodError Rhombohedral{Float64}(1.0, 1.0, 1.0, pi/2, pi/2, pi/3)
    @test Bool(Rhombohedral{Float64}(1.0, 1.0))
    
    @test_throws MethodError Hexagonal{Float64}(2.0, 2.0, 3.0, pi/2, pi/2, pi/3)
    @test Bool(Hexagonal{Float64}(2.0, 3.0))

    @test_throws MethodError Cubic{Float64}(1.0, 1.0, 1.0, pi/2, pi/2, 1/2)
    @test Bool(Cubic{Float64}(1.0))
end

@testset "Peak position calculation" begin
    cl = get_crystal([11.50101, 11.50101, 11.50101, 90.0000, 90.0000, 90.0000], true)
    peak = Peak(2,2,2, 26.853, 100.)
    @test typeof(cl) <: Cubic{Float64}
    @test isapprox(cl(peak), twoθ2q(peak.q), rtol=0.01)

    cl = get_crystal([11.13217, 11.13217, 11.23919, 90.0000, 90.0000,90.0000], true)
    peak = Peak(4, 0, 4, 45.896, 100.)
    @test typeof(cl) <: Tetragonal{Float64}
    @test isapprox(cl(peak), twoθ2q(peak.q), rtol=0.01)

    cl = get_crystal([7.81997, 7.81997, 9.86211, 90.0000, 90.0000, 120.0000], true)
    peak = Peak(1, 1, -2, 29.1958, 100.)
    @test typeof(cl) <: Hexagonal{Float64}
    @test isapprox(cl(peak), twoθ2q(peak.q), rtol=0.01)

    cl = get_crystal([4.26184, 11.11085,10.19209, 90.0000, 90.0000, 90.0000], true)
    peak = Peak(0, 2, 3, 30.851, 100.)
    @test typeof(cl) <: Orthorhombic{Float64}
    @test isapprox(cl(peak), twoθ2q(peak.q), rtol=0.01)

    cl = get_crystal([12.60377, 6.96296, 15.25134, 90.0000, 92.5375, 90.0000], true)
    peak = Peak(2, 2, 1, 30.034, 100.)
    @test typeof(cl) <: Monoclinic{Float64}
    @test isapprox(cl(peak), twoθ2q(peak.q), rtol=0.01)

    cl = get_crystal([6.82618, 7.00339, 7.90972, 106.9500, 89.1806, 101.9832], true)
    peak = Peak(2, -1, 0, 27.3554, 100.)
    @test typeof(cl) <: Triclinic{Float64}
    @test isapprox(cl(peak), twoθ2q(peak.q), rtol=0.01)
    
    peak = Peak(1, 1, -2, 27.9523, 100.)
    @test isapprox(cl(peak), twoθ2q(peak.q), rtol=0.01)
end

@testset "Parameter getters" begin
    cl = Cubic{Float64}(1.)
    @test get_free_params(cl) == [1.]

    cl = Tetragonal{Float64}(1., 3.)
    @test get_free_params(cl) == [1., 3.]

    cl = Hexagonal{Float64}(1., 3.)
    @test get_free_params(cl) == [1., 3.]

    cl = Orthorhombic{Float64}(1., 2., 3.)
    @test get_free_params(cl) == [1., 2., 3.]

    cl = Rhombohedral{Float64}(1., 2.)
    @test get_free_params(cl) == [1., 2.]

    cl = Monoclinic{Float64}(1., 2., 3., 2.5)
    @test get_free_params(cl) == [1., 2., 3., 2.5]

    cl = Triclinic{Float64}(1., 2., 3., 4., 5., 6.)
    @test get_free_params(cl) == [1., 2., 3., 4., 5., 6.]
end