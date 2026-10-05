module TestFullSearch
# Tree search where each node is optimized with full_optimize! (peak heights refit)
using CrystalShift
using CrystalShift: OptimizationSettings, Lazytree, get_phase_ids
using CrystalShift: search!, search_k2n!, TreeSearchSettings, FullOptimizeSettings
using CrystalShift: CrystalPhase, Lorentz, copy_peaks, change_peak_int!
using CrystalShift: Simple, EM, WithUncer
using LinearAlgebra
using Random
using Test

Random.seed!(1)

std_noise = .01
mean_θ = [1., 1., .2]
std_θ = [.5, .5, 1.]

path = "../data/"
f = open(path * "tree_test_sticks.csv", "r")
s = Sys.iswindows() ? split(read(f, String), "#\r\n") : split(read(f, String), "#\n")
close(f)
s[end] == "" && pop!(s)

cs = Vector{CrystalPhase}(undef, size(s))
@. cs = CrystalPhase(String(s), (0.1, ), (Lorentz(), ))
ref_I = [[p.I for p in c.peaks] for c in cs]

x = LinRange(8, 45, 512) |> collect

# Phases 0 and 1 with peak heights scaled away from the reference intensities,
# which a plain optimize! cannot reproduce
function perturbed(c::CrystalPhase)
    c = copy_peaks(c)
    change_peak_int!(c, 0.4 .+ 1.2 .* rand(length(c.peaks)))
    c
end
y = perturbed(cs[1]).(x) + perturbed(cs[2]).(x)
y /= maximum(y)

opt_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 256)
full_opt_stn = FullOptimizeSettings(loop_num=2)

best_node(t) = t[argmin([norm(t[i](x) .- y) for i in eachindex(t)])]
best_residual(t) = minimum(norm(t[i](x) .- y) for i in eachindex(t))

@testset "settings" begin
    d = FullOptimizeSettings()
    @test (d.loop_num, d.mod_peak_num, d.peak_mod_iter, d.analytic_peak_mod) == (8, 32, 32, true)
    @test d.peak_mod_mean == [1.] && d.peak_mod_std == [.5]

    @test isnothing(TreeSearchSettings{Float64}(2, 3, false, false, 5., opt_stn).full_opt_stn)
    @test isnothing(TreeSearchSettings{Float64}().full_opt_stn)
    @test TreeSearchSettings{Float64}(full_opt_stn=full_opt_stn).full_opt_stn === full_opt_stn
    ts = TreeSearchSettings{Float64}(2, 3, false, false, 5., cs[1], opt_stn; full_opt_stn=full_opt_stn)
    @test ts.full_opt_stn === full_opt_stn && ts.default_phase === cs[1]

    uncer_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 256, true, LM, "LS", WithUncer)
    @test_throws ErrorException TreeSearchSettings{Float64}(2, 3, false, false, 5., uncer_stn;
                                                            full_opt_stn=full_opt_stn)
end

@testset "search_k2n! fits perturbed peak heights better than plain optimize!" begin
    plain = search_k2n!(Lazytree(cs, x), x, y, TreeSearchSettings{Float64}(2, 3, false, false, 5., opt_stn))
    full = search_k2n!(Lazytree(cs, x), x, y,
                       TreeSearchSettings{Float64}(2, 3, false, false, 5., opt_stn; full_opt_stn=full_opt_stn))
    @test all(n.is_optimized for n in full)
    @test Set(get_phase_ids(best_node(full))) == Set([0, 1])
    @test best_residual(full) < 0.5 * best_residual(plain)

    # peak heights of the best node moved away from the reference intensities
    bn = best_node(full)
    @test any(any(p.I != ref_I[c.id+1][j] for (j, p) in enumerate(c.peaks)) for c in bn.phase_model.CPs)
    # ... while the phases held by the tree were not modified
    @test all([p.I for p in cs[i].peaks] == ref_I[i] for i in eachindex(cs))
end

@testset "search! levels" begin
    ts = TreeSearchSettings{Float64}(2, 3, false, false, 5., opt_stn; full_opt_stn=full_opt_stn)
    t = search!(Lazytree(cs, x), x, y, ts)
    @test length(t) == 3
    @test all(n.is_optimized for level in t[2:end] for n in level)
    @test all(length(get_phase_ids(n)) == 2 for n in t[3])
    @test Set(get_phase_ids(best_node(reduce(vcat, t[2:end])))) == Set([0, 1])

    # vector k
    ts = TreeSearchSettings{Float64}(3, [3, 2, 1], false, false, 5., opt_stn; full_opt_stn=full_opt_stn)
    t = search!(Lazytree(cs, x), x, y, ts)
    @test length(t) == 4
    @test all(n.is_optimized for level in t[2:end] for n in level)
end

@testset "amorphous root, background, default phase, EM" begin
    # the amorphous root has no phases and falls back to optimize!
    ts = TreeSearchSettings{Float64}(1, 3, true, false, 5., opt_stn; full_opt_stn=full_opt_stn)
    t = search!(Lazytree(cs, x), x, y, ts)
    @test isnothing(t[1][1].phase_model.CPs) && t[1][1].is_optimized
    @test all(n.is_optimized for n in t[2])

    noisy = y .+ 0.1 .* (1 .+ sin.(0.2x))
    ts = TreeSearchSettings{Float64}(2, 2, false, true, 5., opt_stn; full_opt_stn=full_opt_stn)
    t = search_k2n!(Lazytree(cs, x), x, noisy, ts)
    @test all(!isnothing(n.phase_model.background) for n in t)
    res = [norm(t[i](x) .- noisy) for i in eachindex(t)]
    @test Set(get_phase_ids(t[argmin(res)])) == Set([0, 1])

    ts = TreeSearchSettings{Float64}(2, 3, false, false, 5., cs[1], opt_stn; full_opt_stn=full_opt_stn)
    t = search!(Lazytree(cs, x), x, y, ts)
    @test all(0 in get_phase_ids(n) for level in t for n in level)

    em_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 64, true, LM, "LS", EM)
    ts = TreeSearchSettings{Float64}(1, 3, false, false, 5., em_stn; full_opt_stn=FullOptimizeSettings(loop_num=1))
    t = search!(Lazytree(cs, x), x, y, ts)
    @test all(n.is_optimized for n in t[2])
end

println("End of fullsearch.jl test")
end
