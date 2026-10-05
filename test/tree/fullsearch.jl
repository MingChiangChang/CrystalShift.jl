module TestFullSearch
# Tree search whose returned nodes are refit with full_optimize! (peak heights too),
# followed by get_probabilities for the phase assignment
using CrystalShift
using CrystalShift: OptimizationSettings, Lazytree, get_phase_ids
using CrystalShift: search!, search_k2n!, TreeSearchSettings, FullOptimizeSettings
using CrystalShift: CrystalPhase, Lorentz, copy_peaks, change_peak_int!
using CrystalShift: Simple, EM, WithUncer, get_probabilities
using LinearAlgebra
using Test

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
phase(id) = cs[findfirst(c -> c.id == id, cs)]

x = LinRange(8, 45, 512) |> collect

# Test data is built without an RNG so it is the same on every Julia version.
# Peak heights scaled by factors in [0.4, 1.6], which a plain optimize! cannot reproduce
function perturbed(c::CrystalPhase, shift::Real = 0.)
    c = copy_peaks(c)
    change_peak_int!(c, 1 .+ 0.6 .* sin.(2.3 .* (1:length(c.peaks)) .+ shift))
    c
end
# white-looking noise with standard deviation σ (sine hash)
pseudo_noise(n::Int, σ::Real) = σ * sqrt(3) .* (2 .* mod.(sin.(12.9898 .* (1:n)) .* 43758.5453, 1) .- 1)
normalized(y) = y / maximum(y)
reference_data(ids) = normalized(sum(phase(id).(x) for id in ids))
perturbed_data(ids) = normalized(sum(perturbed(phase(id), i - 1.).(x) for (i, id) in enumerate(ids)))

opt_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 256)
full_opt_stn = FullOptimizeSettings(loop_num=2)
settings(f; depth = 2, k = 3, amorphous = false, background = false) =
    TreeSearchSettings{Float64}(depth, k, amorphous, background, 5., opt_stn; full_opt_stn=f)

# phases 3 and 5 with perturbed heights: plain optimize! picks phase 4 instead of 3
y = perturbed_data([3, 5])

residual(n, y) = norm(n(x) .- y)
best_node(t, y) = t[argmin([residual(n, y) for n in t])]
combos(t) = Set(Set(get_phase_ids(n)) for n in t)

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

@testset "search_k2n! refines the searched nodes" begin
    plain = search_k2n!(Lazytree(cs, x), x, y, settings(nothing))
    full = search_k2n!(Lazytree(cs, x), x, y, settings(full_opt_stn))
    # the search itself is unchanged: same phase combinations explored
    @test combos(full) == combos(plain)
    @test all(n.is_optimized for n in full)
    @test Set(get_phase_ids(best_node(plain, y))) != Set([3, 5])
    @test Set(get_phase_ids(best_node(full, y))) == Set([3, 5])
    @test residual(best_node(full, y), y) < 0.5 * residual(best_node(plain, y), y)

    # peak heights of the best node moved away from the reference intensities
    bn = best_node(full, y)
    @test any(any(p.I != ref_I[c.id+1][j] for (j, p) in enumerate(c.peaks)) for c in bn.phase_model.CPs)
    # ... while the phases held by the tree were not modified
    @test all([p.I for p in cs[i].peaks] == ref_I[i] for i in eachindex(cs))
end

@testset "search! levels" begin
    t = search!(Lazytree(cs, x), x, y, settings(full_opt_stn))
    @test length(t) == 3
    @test all(n.is_optimized for level in t[2:end] for n in level)
    @test all(length(get_phase_ids(n)) == 2 for n in t[3])
    plain = search!(Lazytree(cs, x), x, y, settings(nothing))
    @test all(combos(t[l]) == combos(plain[l]) for l in eachindex(t))

    # vector k
    t = search!(Lazytree(cs, x), x, y, settings(full_opt_stn; depth = 3, k = [3, 2, 1]))
    @test length(t) == 4
    @test all(n.is_optimized for level in t[2:end] for n in level)
end

@testset "amorphous root, background, default phase, EM" begin
    # the amorphous root has no phases and is not refit
    t = search!(Lazytree(cs, x), x, y, settings(full_opt_stn; depth = 1, amorphous = true))
    @test isnothing(t[1][1].phase_model.CPs) && t[1][1].is_optimized
    @test all(n.is_optimized for n in t[2])

    noisy = y .+ 0.1 .* (1 .+ sin.(0.2x))
    t = search_k2n!(Lazytree(cs, x), x, noisy, settings(full_opt_stn; k = 2, background = true))
    @test all(!isnothing(n.phase_model.background) for n in t)
    @test Set(get_phase_ids(best_node(t, noisy))) == Set([3, 5])

    ts = TreeSearchSettings{Float64}(2, 3, false, false, 5., cs[1], opt_stn; full_opt_stn=full_opt_stn)
    t = search!(Lazytree(cs, x), x, y, ts)
    @test all(0 in get_phase_ids(n) for level in t for n in level)

    em_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 64, true, LM, "LS", EM)
    ts = TreeSearchSettings{Float64}(1, 3, false, false, 5., em_stn; full_opt_stn=FullOptimizeSettings(loop_num=1))
    t = search!(Lazytree(cs, x), x, y, ts)
    @test all(n.is_optimized for n in t[2])
end

# Full workflow as in paper/calibration.jl: tree search -> drop the root level ->
# get_probabilities -> assign the most probable phase combination.
# Returns the probability of each explored phase combination (duplicates summed).
function phase_probabilities(y, f)
    nodes = reduce(vcat, search!(Lazytree(cs, x), x, y, settings(f))[2:end])
    prob = get_probabilities(nodes, x, y, mean_θ, std_θ)
    @test sum(prob) ≈ 1
    @test all(≥(0), prob)
    p = Dict{Set{Int}, Float64}()
    for (n, pn) in zip(nodes, prob)
        ids = Set(get_phase_ids(n))
        p[ids] = get(p, ids, 0.) + pn
    end
    return p
end
assigned(p) = argmax(p)
p_true(p, ids) = get(p, Set(ids), 0.)

@testset "search + probabilities: phase assignment" begin
    # reference peak heights: refinement keeps the correct assignment
    y_ref = reference_data([1, 3])
    for f in (nothing, full_opt_stn)
        p = phase_probabilities(y_ref, f)
        @test assigned(p) == Set([1, 3])
        @test p_true(p, [1, 3]) > 0.9
    end

    # perturbed peak heights where plain search misassigns (phase 4 for 3)
    for ids in ([3, 5], [3, 7])
        y_pert = perturbed_data(ids)
        p_plain = phase_probabilities(y_pert, nothing)
        p = phase_probabilities(y_pert, full_opt_stn)
        @test assigned(p_plain) != Set(ids)
        @test assigned(p) == Set(ids)
        @test p_true(p, ids) > 0.9
    end

    # perturbed peak heights where plain search is already right
    for ids in ([1, 6], [5, 9])
        p = phase_probabilities(perturbed_data(ids), full_opt_stn)
        @test assigned(p) == Set(ids)
        @test p_true(p, ids) > 0.9
    end

    # one phase with perturbed heights: refinement must not prefer adding a phase
    p = phase_probabilities(perturbed_data([3]), full_opt_stn)
    @test assigned(p) == Set([3])
    @test p_true(p, [3]) > 0.9

    # perturbed peak heights with noise. Not [3, 5]: phases 3 and 4 are both Fd-3m with
    # close lattices, and with noise {4, 5} gets nearly as much probability (0.46 vs 0.52)
    for ids in ([1, 6], [5, 9])
        y_noisy = perturbed_data(ids) .+ pseudo_noise(length(x), 0.02)
        p = phase_probabilities(y_noisy, full_opt_stn)
        @test assigned(p) == Set(ids)
        @test p_true(p, ids) > 0.5
    end
end

println("End of fullsearch.jl test")
end
