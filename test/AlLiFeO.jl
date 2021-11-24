using PhaseMapping: load
using ProgressBars

include("../src/CrystalShift.jl")

std_noise = .01
mean_θ = [1., 1., .2]
std_θ = [.025, 10., 1.]

test_path = "data/AlLiFeO/sticks.csv"
f = open(test_path, "r")
s = split(read(f, String), "#\n") # Windows: #\r\n ...
#s = Vector{CrystalPhase}()
cs = CrystalPhase.(String.(s[1:end-1]))

data, _ = load("AlLiFe", "/Users/mingchiang/Downloads/")


for i in tqdm(eachcol(data.I))
    fit_phases(cs, data.Q, i, std_noise, mean_θ, std_θ)
end