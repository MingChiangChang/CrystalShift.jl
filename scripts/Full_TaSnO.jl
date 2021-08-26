using NPZ
using Plots
using PhaseMapping: xray
using BackgroundSubtraction: mcbl
using LinearAlgebra

# Readin data
include("../src/CrystalShift.jl")
dl = "/Users/mingchiang/Downloads/"
data = npzread(dl * "12_20F16_Ta-Sn-O_integrated.npy")
q = npzread(dl * "12_20F16_Ta-Sn-O_Q.npy")
# x =
# y =
# tpeak =
# dwell =

# CrystalPhas object creation
path = "/Users/mingchiang/Desktop/github/Crystallography_based_shifting/data/"
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

for i in 100:102 # size(data, 1)
    W, H, K = xray(Array(transpose(data[i, :, :])), 4)
    println(size(W), size(H))
    nmf = plot(q[i, :], W)
    display(nmf)
    # plot(H)
    # display(H)

    for j in 1:size(W, 2)
        b = mcbl(W[:, j], q[i,:], 7)
        new = W[:, j] - b
        @. new = max(new, 0)
        println(j)
        @time p = fit_phases(cs, q[i, :], new)
        println(j, p)
        plt = plot(q[i,:], p(q[i, :]))
        plot!(q[i, :], new)
        display(plt)
    end
end
