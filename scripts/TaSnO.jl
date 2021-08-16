using NPZ
using Plots
using PhaseMapping: xray
using BackgroundSubtraction: mcbl

include("../src/CrystalShift.jl")

path = "/Users/mingchiang/Desktop/github/Crystallography_based_shifting/data/"

Q = npzread(path * "Test_Ta2O5_Q.npy")
data = npzread(path * "Test_Ta2O5_XRD_map.npy")

# CrystalPhas object creation
phase_path = path * "Ta-Sn-O/sticks.csv"
f = open(phase_path, "r")
s = split(read(f, String), "#\n") # Windows: #\r\n ...
if s[end] == ""
    pop!(s)
end
cs = Vector{CrystalPhase}(undef, size(s))
for i in eachindex(s)
    println(i)
    cs[i] = CrystalPhase(String(s[i]))
end

x = collect(12:.01:40)
# plot(x, cs[6](x), label="ICSD", linewidth=1.4)
# plot!(x, reconstruct!(cs[6], [6.24, 3.62, 7.7, 1., .1], x), label="Modified", linewidth=1.4)
# xlabel!("Q (1/nm)")
# ylabel!("a. u.")
# title!("Ta2O5")
# savefig("ta2o5.png")
# heatmap(data)
W, H, K = xray(data, 4)
# plot(H')

for i in eachcol(W)
    plt = plot()
    b = mcbl(i, Q, 7)
    new = i-b
    @. new = max(new, 0)
    plot!(Q, new)

    phase = optimize!(cs[6], Q, new; regularization=false)
    plot!(Q, phase(Q))
    display(plt)
end
