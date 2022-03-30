using ProgressBars

using CrystalShift
using CrystalShift: get_free_params, extend_priors, Lorentz, evaluate_residual!, PseudoVoigt
using CrystalShift: Gauss, FixedPseudoVoigt, PeakModCP, optimize!, evaluate!, full_optimize!
using PhaseMapping: load
using Plots
using LinearAlgebra
using BenchmarkTools

std_noise = 1e-2
# mean_θ = [1., 1., .1, 1.] # Set to favor true solution
# std_θ = [.2, .5, .2, 1.]

mean_θ = [1., 1., .1] # Set to favor true solution
std_θ = [0.2, .5, 1.]

peak_mean_θ = [1.]
peak_std_θ = [1.]

method = LM
objective = "LS"
improvement = 0.1

test_path = "/Users/ming/Downloads/AlLiFeO/sticks.csv"
# test_path = "/Users/ming/Downloads/cif/sticks.csv"
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n") # Windows: #\r\n ...
else
    s = split(read(f, String), "#\n")
end

if s[end] == ""
    pop!(s)
end

cs = Vector{CrystalPhase}(undef, size(s))
# cs = @. PeakModCP(String(s), (0.1, ), (FixedPseudoVoigt(0.01), ))
cs = @. CrystalPhase(String(s), (0.1, ), (FixedPseudoVoigt(0.01), ))
# println("$(size(cs, 1)) phase objects created!")
max_num_phases = 3
data, _ = load("AlLiFe", "/Users/ming/Downloads/")
x = data.Q
x = x[1:400]

y = data.I[1:400, 1]
y /= maximum(y)

# c = optimize()
@time c = full_optimize!([cs[6], ], x, y, std_noise, mean_θ, std_θ;
                  objective = objective, method = method, maxiter = 128,
                  regularization = true, verbose = false)

# pm_cp = PeakModCP.(c)


# @time cp = optimize!(pm_cp, x, y, std_noise, peak_mean_θ, peak_std_θ;
#                   objective = objective, method = bfgs, maxiter = 16,
#                   regularization = true, verbose = false)