# Phase assignment through the tree search with method = LM vs dogleg.
#
#   julia --project=benchmark -t 4 benchmark/tree.jl
#
# 7 phase pairs from data/tree_test_sticks.csv × perturbed / reference peak heights ×
# with / without noise (28 patterns), each run through search! (depth 2, k 3), optionally
# with FullOptimizeSettings refinement, then get_probabilities; a pattern counts as correct
# when the most probable phase combination is the true one. Data is built without an RNG.
using Pkg
Pkg.develop(path = joinpath(@__DIR__, ".."); io = devnull)
Pkg.instantiate(; io = devnull)

using CrystalShift
using CrystalShift: OptimizationSettings, Lazytree, get_phase_ids, search!, TreeSearchSettings, FullOptimizeSettings
using CrystalShift: CrystalPhase, Lorentz, copy_peaks, change_peak_int!, get_probabilities, LM, dogleg
using Printf

std_noise, mean_θ, std_θ = .01, [1., 1., .2], [.5, .5, 1.]
blocks = filter(!isempty ∘ strip, split(read(joinpath(@__DIR__, "..", "data", "tree_test_sticks.csv"), String), "#\n"))
cs = CrystalPhase.(String.(blocks), (0.1,), (Lorentz(),))
x = collect(LinRange(8, 45, 512))
phase(id) = cs[findfirst(c -> c.id == id, cs)]
function perturbed(c, shift)  # peak heights scaled by factors in [0.4, 1.6]
    c = copy_peaks(c); change_peak_int!(c, 1 .+ 0.6 .* sin.(2.3 .* (1:length(c.peaks)) .+ shift)); c
end
pseudo_noise(n, σ) = σ * sqrt(3) .* (2 .* mod.(sin.(12.9898 .* (1:n)) .* 43758.5453, 1) .- 1)

function assign(y, method, refine)
    opt = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 256, true, method)
    ts = TreeSearchSettings{Float64}(2, 3, false, false, 5., opt; full_opt_stn = refine ? FullOptimizeSettings() : nothing)
    t = @elapsed nodes = reduce(vcat, search!(Lazytree(cs, x), x, y, ts)[2:end])
    p = get_probabilities(nodes, x, y, mean_θ, std_θ)
    Set(get_phase_ids(nodes[argmax(p)])), t
end

println("Julia ", VERSION, ", threads ", Threads.nthreads(), ", dogleg lattice bound ±", 100CrystalShift.DOGLEG_MAX_STRAIN, "%")
for method in (LM, dogleg); assign(phase(1).(x) .+ phase(3).(x), method, true); end  # compile
for refine in (false, true), method in (LM, dogleg)
    correct, time, n = 0, 0., 0
    for ids in ([1,3],[3,5],[1,6],[5,9],[6,12],[1,4],[3,7]), pert in (true, false), σ in (0., .02)
        y = pert ? sum(perturbed(phase(id), i - 1.).(x) for (i, id) in enumerate(ids)) : sum(phase(id).(x) for id in ids)
        y = y ./ maximum(y) .+ pseudo_noise(length(x), σ)
        a, t = assign(y, method, refine)
        correct += a == Set(ids); time += t; n += 1
    end
    @printf("%-7s refine=%-5s correct %2d/%d   total search time %.1f s\n", method, refine, correct, n, time)
end
