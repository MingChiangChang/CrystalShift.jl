# Run the benchmark suite and save results.
#
#   julia --project=benchmark -t 4 benchmark/run.jl [name] [--compare baseline] [--filter substring]
#
# Results go to benchmark/results/<name>.json (default name: short git hash, with
# "-dirty" if the tree has uncommitted changes). With --compare, each benchmark is
# judged against benchmark/results/<baseline>.json (median time, 5% tolerance).
using Pkg
Pkg.develop(path = joinpath(@__DIR__, ".."); io = devnull)
Pkg.instantiate(; io = devnull)

using BenchmarkTools
using Printf

function parse_args(args)
    name, baseline, filt = nothing, nothing, nothing
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--compare"
            baseline = args[i += 1]
        elseif a == "--filter"
            filt = args[i += 1]
        else
            name = a
        end
        i += 1
    end
    if name === nothing
        cd(@__DIR__) do
            hash = readchomp(`git rev-parse --short HEAD`)
            dirty = !isempty(readchomp(`git status --porcelain -- ../src`))
            name = dirty ? "$hash-dirty" : hash
        end
    end
    name, baseline, filt
end

name, baseline, filt = parse_args(ARGS)
resultdir = joinpath(@__DIR__, "results")
mkpath(resultdir)

include(joinpath(@__DIR__, "benchmarks.jl"))

suite = SUITE
if filt !== nothing
    keep = BenchmarkGroup()
    for (k, b) in BenchmarkTools.leaves(SUITE)
        occursin(filt, join(k, "/")) && (keep[k] = b)
    end
    suite = keep
end

println("Threads: ", Threads.nthreads(), "   Julia: ", VERSION)
println("Correctness checks:")
for (k, v) in sort(collect(checks()); by = first)
    println("  ", rpad(k, 36), v)
end

println("\nTuning and running…")
tune!(suite)
results = run(suite; verbose = false, seconds = 10)
med = median(results)

println()
for (k, t) in sort(collect(BenchmarkTools.leaves(med)); by = first)
    @printf("  %-45s %12s  %10s  %8d allocs\n", join(k, " / "), BenchmarkTools.prettytime(time(t)),
            BenchmarkTools.prettymemory(memory(t)), allocs(t))
end

outfile = joinpath(resultdir, "$name.json")
BenchmarkTools.save(outfile, results)
println("\nSaved to ", relpath(outfile, pwd()))

if baseline !== nothing
    base = BenchmarkTools.load(joinpath(resultdir, "$baseline.json"))[1]
    println("\nComparison against $baseline (median time, ratio < 1 is faster):")
    for (k, t) in sort(collect(BenchmarkTools.leaves(med)); by = first)
        b = try base[k] catch; continue end
        j = judge(t, median(b); time_tolerance = 0.05)
        @printf("  %-45s %6.2fx  %s\n", join(k, " / "), ratio(j).time, time(j))
    end
end
