module Testutil

using CrystalShift: parse_cond, get_weighted_center, kl
using CrystalShift: negative_log_poisson, negative_log_poisson_λ, poisson
using CrystalShift: Gauss
using Test

g = Gauss()
x = collect(-5:.01:5)
y = [g.(x) for _ in 1:5]
y = mapreduce(permutedims, vcat, y)

@testset "util" begin
    @test parse_cond("chess_x_y_700_1000", Float64) == [700., 1000.]
    @test isapprox(get_weighted_center(y), 500, atol=2)
    @test kl(y, y) == 0.
    @test kl(1, 3) ≈ -1.0986122886681098
    @test poisson(1, 1) ≈ 0.36787944117144233
    @test negative_log_poisson(1, 1) ≈ -log(0.36787944117144233)
    @test negative_log_poisson(1.2, 1) ≈ -log(0.36143305429464256)
    @test negative_log_poisson_λ(1, 1) ≈ -log(0.36787944117144233)
end



end # module