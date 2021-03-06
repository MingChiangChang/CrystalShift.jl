using NPZ
using Plots
using PhaseMapping: xray
using BackgroundSubtraction: mcbl
using LinearAlgebra
using JSON
using TimerOutputs
using ProgressBars
using Profile
using CSV
using DataFrames

# const to = TimerOutput()

# Constant Declaration
const RANK = 4
const THRESH = 0.3

# Readin data
include("../src/CrystalShift.jl")
# dl = "/Users/mingchiang/Downloads/"
dl = "/Users/r2121/Downloads/"
# data = npzread(dl * "12_20F16_Ta-Sn-O_integrated.npy")
# q = npzread(dl * "12_20F16_Ta-Sn-O_Q.npy")
data = npzread(dl * "TaSnO_data.npy")
q = npzread(dl * "TaSnO_Q.npy")

f = open(dl * "TaSnO_conds.json", "r")
cond = JSON.parse(f)
conds = parse_cond.(cond, Float64) # [x, y, tpeak, dwell]

# Load composition
df = DataFrame(CSV.File(dl * "59778_TaSn_20F16_DwellTpeak.csv"))
ta = df[!, "Ta.nmol_offgrid"]
sn = df[!, "Sn.nmol_offgrid"]
cation = zero(ta)
@. cation = ta/(ta+sn)

# CrystalPhas object creation
path = "data/"
phase_path = path * "Ta-Sn-O/sticks.csv"
f = open(phase_path, "r")
s = split(read(f, String), "#\r\n") # Windows: #\r\n ...

if s[end] == ""
    pop!(s)
end

cs = Vector{CrystalPhase}(undef, size(s))
for i in eachindex(s)
    cs[i] = CrystalPhase(String(s[i]))
end
println("$(size(cs)) phase objects created!")

wafer_result = Vector{StripeResult}(undef, size(data)[1])

for i in tqdm(1:size(data, 1)) # size(data, 1)
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

    stripe = Vector{PhaseResult}(undef, 3)
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
        p = fit_phases(cs, q[i, :], new)
        # plt = plot(q[i,:], p(q[i, :]))
        # plot!(q[i, :], new)
        # display(plt)

        stripe[j] =  PhaseResult(p.cl, p.name, BH[j, :], new, isCenter[j])
    end
    subdf = df[(df.xcenter .== conds[i][1]) .& (df.ycenter .== conds[i][2]), :]
    local ta = subdf[!, "Ta.nmol_offgrid"]
    local sn = subdf[!, "Sn.nmol_offgrid"]
    local cation = ta./(ta.+sn)
    wafer_result[i] =  StripeResult(stripe, cation[1], conds[i]...)
end
# Does data make sense

# Store data
open("data/TaSnO.json", "w") do f
    JSON.print(f, wafer_result)
end

# Plotting
