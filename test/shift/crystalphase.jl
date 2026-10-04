# include("../src/CrystalShift.jl")
module Testcrystalshift
using Test
using CrystalShift
using CrystalShift: evaluate!, get_param_nums, get_eight_params, get_free_lattice_params
using CrystalShift: collect_crystals, Lorentz, PseudoVoigt, get_moles, Gauss
using CrystalShift: get_fraction, get_free_params, evaluate_residual!, get_intrinsic_profile_type
using NPZ
using LinearAlgebra

path = "../data/"
test_path = path * "Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")
sol = npzread(path * "crystalphase_test_sol.npy")
multi_sol = npzread(path * "crystalphase_multi_sol.npy")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]), (0.1,), (Lorentz(),))
f = open(test_path, "r")
cs2 = CrystalPhase(f)
gauss_cs = CrystalPhase(String(s[1]), 0.1, Gauss())
x = collect(8:.1:60)
y = zero(x)

show(cs[1])

@testset "Helper functions" begin
    @test Bool(cs[1]) == true
    @test Bool(cs) == true
    @test get_param_nums(cs[1]) == 6
end

@testset "CrystalPhase object creation" begin
    @test Bool(CrystalPhase(String(s[1])))
    @test Bool(CrystalPhase.(String.(s[1:end-1])))
end

@testset "Getters" begin
    @test get_eight_params(cs[1]) ≈ [17.11299, 4.872, 5.548, deg2rad(90.0), deg2rad(90.59), deg2rad(90.0), 1.0, 0.1]
    @test get_eight_params(cs[1], [1., 2., 3., 4., 5., 6.]) == [1., 2., 3., pi/2, 4., pi/2, 5., 6.]
    @test get_eight_params(cs[2], [1., 2., 3.]) == [1., 1., 1., pi/2, pi/2, pi/2, 2., 3.]
    @test get_eight_params(cs[6], [1., 2., 3., 4., 5.]) == [1., 2., 3., pi/2, pi/2, pi/2, 4., 5.]
    @test get_eight_params(cs[10], [1., 2., 3., 4.]) == [1., 1., 2., pi/2, pi/2, pi/2, 3., 4.]

    @test get_free_lattice_params(cs[1]) ≈ [17.11299, 4.872, 5.548, deg2rad(90.59)]
    @test get_free_lattice_params(cs[1], [1., 2., 3., 4., 5., 6.]) == [1., 2., 3., 5.]
    @test get_free_lattice_params(cs[2], [1., 2., 3., 4., 5., 6.]) == [1.]
    @test get_free_lattice_params(cs[6], [1., 2., 3., 4., 5., 6.]) == [1., 2., 3.]
    @test get_free_lattice_params(cs[10], [1., 2., 3., 4., 5., 6.]) == [1., 3.]

    @test collect_crystals(cs[1:3]) == [cs[1].cl, cs[2].cl, cs[3].cl]
    @test show(cs[1])
    new_cs2 = CrystalPhase(cs[2], [1., 2., 3.,])
    new_cs6 = CrystalPhase(cs[6], [1., 2., 3., 4., 5.])
    @test get_moles(new_cs2) ≈ 0.06429554942685066
    @test get_moles(new_cs6) ≈ 0.59129357407787
    @test get_fraction([new_cs2, new_cs6]) ≈ [get_moles(new_cs2), get_moles(new_cs6)]./(get_moles(new_cs2)+get_moles(new_cs6))
    @test get_intrinsic_profile_type(typeof(gauss_cs.profile)) == Gauss
end

# TODO: Generate new test sets
# Test construction with 5 phases from different crystal family
@testset "Single phase evaluation" begin
    pn = [1,2,5,6,10]
    for i in eachindex(pn)
        e = zero(x)
        sol[i,:] ./= maximum(sol[i,:])
        evaluate!(e, cs[pn[i]], x)
        e ./= maximum(e)
        # println(maximum(sol[i,:] .-e))
        @test isapprox(sol[i,:], e, atol=0.05)
    end
end

@testset "Multiphase evaluation" begin
    pn = [[1,2], [1,5], [2,5], [5,6], [6,10]]
    for i in eachindex(pn)
        e = zero(x)
        multi_sol[i,:] ./= maximum(multi_sol[i,:])
        evaluate!(e, cs[pn[i]], x)
        e ./= maximum(e)
        @test isapprox(multi_sol[i,:], e, atol=0.05)
    end
end

t = zero(x)
evaluate!(t, cs[1], x)
# println(get_free_params(cs[1]))
evaluate_residual!(cs[1], get_free_params(cs[1]), x, t)
@test norm(t) < 10^-10

t = zero(x)
evaluate!(t, cs[1], get_free_params(cs[1]), x)
# println(get_free_params(cs[1]))
evaluate_residual!(cs[1], get_free_params(cs[1]), x, t)
@test norm(t) < 10^-10

t = zero(x)
evaluate!(t, cs[1:2], get_free_params(cs[1:2]), x)
# println(get_free_params(cs[1:2]))
evaluate_residual!(cs[1:2], get_free_params(cs[1:2]), x, t)
@test norm(t) < 10^-10

# cs = CrystalPhase.(String.(s[1:end-1]), (0.1,), (PseudoVoigt(0.5),))
# x = collect(8:.1:60)
# y = zero(x)

# @testset "Helper functions Voigt" begin
#     @test Bool(cs[1]) == true
#     @test Bool(cs) == true
#     @test get_param_nums(cs[1]) == 7
# end

# @testset "Getters Voigt" begin
#     @test get_eight_params(cs[1]) ≈ [17.11299, 4.872, 5.548, deg2rad(90.0), deg2rad(90.59), deg2rad(90.0), 1.0, 0.1, 0.5]
#     @test get_eight_params(cs[1], [1., 2., 3., 4., 5., 6.]) == [1., 2., 3., pi/2, 4., pi/2, 5., 6., 0.5]
#     @test get_eight_params(cs[2], [1., 2., 3.]) == [1., 1., 1., pi/2, pi/2, pi/2, 2., 3., 0.5]
#     @test get_eight_params(cs[6], [1., 2., 3., 4., 5.]) == [1., 2., 3., pi/2, pi/2, pi/2, 4., 5., 0.5]
#     @test get_eight_params(cs[10], [1., 2., 3., 4.]) == [1., 1., 2., pi/2, pi/2, pi/2, 3., 4., 0.5]

#     @test get_free_lattice_params(cs[1]) ≈ [17.11299, 4.872, 5.548, deg2rad(90.59)]
#     @test get_free_lattice_params(cs[1], [1., 2., 3., 4., 5., 6.]) == [1., 2., 3., 5.]
#     @test get_free_lattice_params(cs[2], [1., 2., 3., 4., 5., 6.]) == [1.]
#     @test get_free_lattice_params(cs[6], [1., 2., 3., 4., 5., 6.]) == [1., 2., 3.]
#     @test get_free_lattice_params(cs[10], [1., 2., 3., 4., 5., 6.]) == [1., 3.]

#     @test collect_crystals(cs[1:3]) == [cs[1].cl, cs[2].cl, cs[3].cl]
# end

# # Test construction with 5 phases from different crystal family
# @testset "Single phase evaluation Voigt" begin
#     pn = [1,2,5,6,10]
#     for i in eachindex(pn)
#         e = zero(x)
#         sol[i,:] ./= maximum(sol[i,:])
#         evaluate!(e, cs[pn[i]], x)
#         e ./= maximum(e)
#         @test sol[i,:] ≈ e
#     end
# end

# @testset "Multiphase evaluation Voigt" begin
#     pn = [[1,2], [1,5], [2,5], [5,6], [6,10]]
#     for i in eachindex(pn)
#         e = zero(x)
#         multi_sol[i,:] ./= maximum(multi_sol[i,:])
#         evaluate!(e, cs[pn[i]], x)
#         e ./= maximum(e)
#         @test multi_sol[i,:] ≈ e
#     end
# end


end # modules