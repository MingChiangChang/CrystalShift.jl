module TestPeakProfile
using Test
using PhaseMapping: PeakProfile, Lorentz, Gauss, PseudoVoigt

@testset "PeakProfile" begin
    l = Lorentz()
    g = Gauss()
    v = PseudoVoigt(.5)

    x = randn()
    μ = randn()
    σ = exp(randn())

    for p in (l, g, v)
        @test p isa PeakProfile
        @test p(x, μ, σ) ≈ p((x-μ)/σ)
    end
    @test 2v(x) ≈ g(x) + l(x) 
end

end # TestPeakProfile