module Testtree
using CrystalShift
using CrystalShift: search!, bft, pos_res_thresholding, get_level, get_phase_ids
# using DelimitedFiles
using LinearAlgebra
using Test

using CrystalShift: CrystalPhase, optimize!, Lorentz, PseudoVoigt

std_noise = .01
mean_θ = [1., 1., 1.]
std_θ = [.5, 0.5, 1.]

# CrystalPhas object creation
path = "../data/"
phase_path = path * "tree_test_sticks.csv"
f = open(phase_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n") # Windows: #\r\n ...
else
    s = split(read(f, String), "#\n")
end

if s[end] == ""
    pop!(s)
end

cs = Vector{CrystalPhase}(undef, size(s))
@. cs = CrystalPhase(String(s), (0.1, ), (Lorentz(), ))

x = LinRange(8, 45, 1024)
y = cs[1].(x)+cs[2].(x)
y /= max(y...)

# TODO: add unit test for smaller function

a = Tree(cs, 2)

traversal = bft(a)
@testset "bft test" begin
    @test all(n->get_level(n) == 1, traversal[1:15])
    @test all(n->get_level(n) == 2, traversal[16:120])
end

# TODO: find bound error
@time res = search!(a, bft, x, y, std_noise,
                    mean_θ, std_θ, pos_res_thresholding,
                    maxiter = 64, regularization = true, tol = 1.)
residual = Float64[]

ind = argmin([norm(i(x)-y) for i in res])

@testset "Basic tree search with trimming" begin
    @test Set(get_phase_ids(res[ind])) == Set([0, 1])
end
println("End of tree.jl test")
end # module

