using NPZ
using Plots
using PhaseMapping: xray
using BackgroundSubtraction: mcbl

include("../src/CrystalShift.jl")
#using CrystalShift: CrystalPhase, optimize!

path = "/Users/mingchiang/Desktop/github/CrystalShift/data/"
#path = "/Users/r2121/Desktop/Code/Crystallography_based_shifting/data/"
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
    cs[i] = CrystalPhase(String(s[i]))
end

x = collect(12:.01:40)
h = heatmap(data)
display(h)

W, H, K = xray(data, 4)
nmf = plot(Q, W)
display(nmf)

for i in eachcol(W)
    plt = plot()
    b = mcbl(i, Q, 7)
    new = i-b
    @. new = max(new, 0)
    plot!(Q, new)

    #@time phase = optimize!(cs[6], Q, new; regularization=true)
    @time phase = fit_phases(cs, Q, new; regularization=true)
    rec = phase.(Q)
    plot!(Q, rec)
    xlabel!("Q (1/nm)")
    ylabel!("a. u.")
    display(plt)
    # println(norm(rec.-new)/max.(rec...)) why is this so slow?
    println(phase.cl)
end
