module TestPeakProfile
using Test
using CrystalShift: PeakProfile, Lorentz, Gauss, PseudoVoigt
using CrystalShift: get_param_nums, get_free_params, sigmoid

@testset "PeakProfile" begin
    l = Lorentz()
    g = Gauss()
    v = PseudoVoigt{Float64}(.5)

    x = randn()
    μ = randn()
    σ = exp(randn())

    @test get_param_nums(l) == 0
    @test get_param_nums(g) == 0
    @test get_param_nums(v) == 1
    @test get_free_params(l) == []
    @test get_free_params(g) == []
    @test get_free_params(v) == [sigmoid(.5)]

    for p in (l, g, v)
        @test p isa PeakProfile
        @test p(x, μ, σ) ≈ p((x-μ)/σ)
    end
    # @test 2v(x) ≈ g(x) + l(x) 
end

end # TestPeakProfile