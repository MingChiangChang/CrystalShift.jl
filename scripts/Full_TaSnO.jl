using NPZ
using Plots
using PhaseMapping: xray
using BackgroundSubtraction: mcbl
using LinearAlgebra
using JSON
using TimerOutputs

const to = TimerOutput()

# Constant Declaration
RANK = 4

# Readin data
include("../src/CrystalShift.jl")
dl = "/home/mingchiang/Downloads/"
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

for i in 100:101 # size(data, 1)
    W, H, K = xray(Array(transpose(data[i, :, :])), RANK)
    println(size(W), size(H))
    nmf = plot(q[i, :], W)
    display(nmf)
    wanted = collect(1:RANK)
    deleteat!(wanted, argmax(H[:,1]))
    BW = W[:, wanted]
    BH = H[wanted, :]

    for j in 1:size(BW, 2)
        b = mcbl(BW[:, j], q[i,:], 7)
        new = BW[:, j] - b
        @. new = max(new, 0)
        println(j)
        @timeit to "fitting $(i)th stripe $(j)th pattern" p = fit_phases(cs, q[i, :], new)
        println(j, p)
        plt = plot(q[i,:], p(q[i, :]))
        plot!(q[i, :], new)
        display(plt)
    end
end

show(to)