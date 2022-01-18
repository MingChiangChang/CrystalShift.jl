# include("../src/CrystalShift.jl")
# TODO: if this should be in the testing suite, add tests and remove plots and @time macros
module Testcrystalshift
using Test
using CrystalShift
using CrystalShift: evaluate!, get_param_nums, get_eight_params, get_free_lattice_params
using CrystalShift: collect_crystals
using NPZ

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

cs = CrystalPhase.(String.(s[1:end-1]))
x = collect(8:.1:60)

x = collect(8:.1:60)
y = zero(x)

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
end

# Test construction with 5 phases from different crystal family
@testset "Single phase evaluation" begin
    pn = [1,2,5,6,10]
    for i in eachindex(pn)
        e = zero(x)
        @test sol[i,:] ≈ evaluate!(e, cs[pn[i]], x)
    end
end

@testset "Multiphase evaluation" begin
    pn = [[1,2], [1,5], [2,5], [5,6], [6,10]]
    for i in eachindex(pn)
        e = zero(x)
        @test multi_sol[i,:] ≈ evaluate!(e, cs[pn[i]], x)
    end
end

# @test "Single phase evaluation with modified parameters" begin
#     # take a phase then multiply its lattice params by 1.002

# end

# function synthesize_data(cp::CrystalPhase, x::AbstractVector)
#     params = get_free_params(cp)
#     interval_size = 0.025
#     scaling = (interval_size.*rand(size(params, 1),) .- interval_size/2) .+ 1
#     @. params = params*scaling
#     params = [params..., 1., 0.2]
#     r = reconstruct!(cp, params, x)
#     normalization = maximum(r)
#     params[end-1] /= normalization
#     if verbose
#         println("synthesize_data: ", params)
#     end

#     return r/normalization, params
# end

end # modules