module TestLazyTree
using CrystalShift
using CrystalShift: OptimizationSettings
using CrystalShift: Lazytree, get_phase_ids
using CrystalShift: add_phase, create_child_nodes, attach_child_nodes!, expand!
using CrystalShift: search!, search_k2n!, TreeSearchSettings
# using DelimitedFiles
using LinearAlgebra

using CrystalShift: CrystalPhase, optimize!, Lorentz, PseudoVoigt

using Test

std_noise = .01
mean_θ = [1., 1., .2]
std_θ = [.5, .5, 1.]

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
@. cs = CrystalPhase(String(s), (0.1, ), (Lorentz(), )) # For ease of testing fast

x = LinRange(8, 45, 512) |> collect
y = cs[1].(x)+cs[2].(x)
y /= max(y...)
#collect(x)
LT = Lazytree(cs, x)

# test_pm1 = add_phase(LT.nodes[1].phase_model, cs[1], x, 5)
# test_pm2 = add_phase(test_pm1, cs[2], x, 5)


# @test test_pm1.CPs[1] == cs[1]
# @test test_pm2.CPs == cs[1:2]

# child_nodes = create_child_nodes(LT, LT.nodes[1], 1)
# @test get_phase_ids.(child_nodes) == [[i] for i in 0:14]
# attach_child_nodes!(LT.nodes[1], child_nodes)
# @test get_phase_ids.(LT.nodes[1].child_node) == [[i] for i in 0:14]
# push!(LT.nodes, child_nodes...)
# t = expand!(LT, LT.nodes[1].child_node[2])

LT = Lazytree(cs,  x)

noise_intensity = 0.1
noise = noise_intensity.*(1 .+ sin.(0.2x))
# @. y += noise
opt_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 256)
ts_stn = TreeSearchSettings{Float64}(2, 3, false, false, 5., opt_stn)

t = search!(LT, x, y, ts_stn)
ts_stn = TreeSearchSettings{Float64}(2, 3, true, false, 5., opt_stn)
t = search!(LT, x, y, ts_stn)
#t = search!(LT, x, y, 2, 10, std_noise, mean_θ, std_θ,
#                  maxiter=64, regularization=true)
LT = Lazytree(cs, x)
ts_stn = TreeSearchSettings{Float64}(2, 3,  false, false, 5., opt_stn)
@time t = search_k2n!(LT, x, y, ts_stn) # 2.3 secs without background modeling
res = [norm(t[i](x).-y) for i in eachindex(t)]
ind = argmin(res)
@test Set(get_phase_ids(t[ind])) == Set([0, 1])

LT = Lazytree(cs, x)
ts_stn = TreeSearchSettings{Float64}(3, [3, 2, 1],  false, false, 5., opt_stn)
@time t = search!(LT, x, y, ts_stn) # 2.3 secs without background modeling

@. y += noise
LT = Lazytree(cs, x)
ts_stn = TreeSearchSettings{Float64}(2, 3, false, true, 5., opt_stn)
@time t = search_k2n!(LT, x, y, ts_stn) # 28 secs with background modeling
# t = reduce(vcat, t)
res = [norm(t[i](x).-y) for i in eachindex(t)]
ind = argmin(res)
# using Plots
# plt = plot(x, y)
# plot!(x, evaluate!(zero(x), t[ind].phase_model, x))
# display(plt)
@test Set(get_phase_ids(t[ind])) == Set([0, 1])

LT = Lazytree(cs, x)
ts_stn = TreeSearchSettings{Float64}(3, [3, 2, 1],  false, false, 5., cs[1], opt_stn)
@time t = search!(LT, x, y, ts_stn) 

println("End of lazytree.jl test")
end