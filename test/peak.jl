# test peak construction
using Test

include("../src/peak.jl")

Base.Bool(p::Peak) = true

@testset "Peak creation" begin
    @test Bool(Peak(1, 0, 0, 10., 10.))
    @test_throws InexactError Peak(2.5, 2.5, 2, 3, 3)
    @test Bool(Peak("1,2,3,4,5"))
end
