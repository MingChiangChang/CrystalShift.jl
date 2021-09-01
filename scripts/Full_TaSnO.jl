using NPZ
using Plots
using PhaseMapping: xray
using BackgroundSubtraction: mcbl
using LinearAlgebra
using JSON
using TimerOutputs
using ProgressBars
using Profile

const to = TimerOutput()

# Constant Declaration
RANK = 4
THRESH = 0.3

# Readin data
include("../src/CrystalShift.jl")
dl = "/Users/mingchiang/Downloads/"
data = npzread(dl * "12_20F16_Ta-Sn-O_integrated.npy")
q = npzread(dl * "12_20F16_Ta-Sn-O_Q.npy")

f = open(dl * "12_20F16_Ta-Sn-O_cond.json", "r")
cond = JSON.parse(f)
conds = parse_cond.(cond, Float64) # [x, y, tpeak, dwell]

# CrystalPhas object creation
path = "data/"
phase_path = path * "Ta-Sn-O/sticks.csv"
f = open(phase_path, "r")
s = split(read(f, String), "#\n") # Windows: #\r\n ...

if s[end] == ""
    pop!(s)
end

cs = Vector{CrystalPhase}(undef, size(s))
for i in eachindex(s)
    cs[i] = CrystalPhase(String(s[i]))
end
println("$(size(cs)) phase objects created!")

wafer_result = Vector{StripeResult}()

for i in tqdm(100:102) # size(data, 1)
    # TODO Pre-screening of the heatmap
    # TODO Try t-SNE or UMAP on the data?
    W, H, K = xray(Array(transpose(data[i, :, :])), RANK)
    # println(size(W), size(H))
    # nmf = plot(q[i, :], W)
    # display(nmf)
    wanted = collect(1:RANK)
    deleteat!(wanted, argmax(H[:,1]))
    BW = W[:, wanted]
    BH = H[wanted, :]

    stripe = Vector{PhaseResult}()
    center = get_weighted_center(BH)
    isCenter = BitArray(undef, RANK-1)
    for k in 1:3
        isCenter[k] = BH[k, center[k]] > THRESH
    end

    for j in 1:size(BW, 2)
        # TODO Pre-screening of the spectrum
        b = mcbl(BW[:, j], q[i,:], 7)
        new = BW[:, j] - b
        @. new = max(new, 0)
        @timeit to "fitting $(i)th stripe $(j)th pattern" p = fit_phases(cs, q[i, :], new)
        # plt = plot(q[i,:], p(q[i, :]))
        # plot!(q[i, :], new)
        # display(plt)

        push!(stripe, PhaseResult(p, BH[j, :], new, isCenter[j]))
    end
    push!(wafer_result, StripeResult(stripe, [1,1], conds[i]...))
end
# Does data make sense
show(to)
