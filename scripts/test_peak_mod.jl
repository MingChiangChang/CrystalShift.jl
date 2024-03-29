using ProgressBars

using CrystalShift
using CrystalShift: get_free_params, extend_priors, Lorentz, evaluate_residual!, PseudoVoigt
using CrystalShift: Gauss, FixedPseudoVoigt, optimize!, evaluate!, full_optimize!
using CrystalShift: PhaseModel, PeakModCP, get_PeakModCP
# using PhaseMapping: load
using Plots
using LinearAlgebra
using BenchmarkTools
using CovarianceFunctions: EQ
using NPZ

using Base.Threads

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
# data, _ = load("AlLiFe", "/Users/ming/Downloads/")
t = npzread("alfeli_noise=5e-2.npy")
# x = data.Q
# x = x[1:400]
x = collect(15:.1:79.9)

y = t[1, :]
y /= maximum(y)

# for i in eachindex(cs)
# c = optimize()
bg = BackgroundModel(x, EQ(), 5, rank_tol=1e-2)
# @time c = full_optimize!(PhaseModel(cs[17], bg), x, y, std_noise, mean_θ, std_θ;
#                 objective = objective, method = method, maxiter = 128,
#                 regularization = true, verbose = true, tol=1E-5)
@time c = full_optimize!(PhaseModel(cs[17], bg), x, y, std_noise, mean_θ, std_θ;
                            mod_peak_num = 200, loop_num=32,
                            method=LM, peak_shift_iter=256, peak_mod_iter=16,
                            objective="LS", regularization=true, peak_mod_mean=[1., 1.] , peak_mod_std=[.5, .05]
)
plot(x, y)
plt = plot(x, y)
plot!(x, evaluate!(zero(x), c, x))
display(plt)
# end
# @time pm_cp = PeakModCP.(c.CPs, (x, ), (32, ))
# @time pm_cp = get_PeakModCP(c, x, 32)

# @time cp = optimize!(pm_cp, x, y, std_noise, peak_mean_θ, peak_std_θ;
#                   objective = objective, method = bfgs, maxiter = 16,
#                   regularization = true, verbose = false)