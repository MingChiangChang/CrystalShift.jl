module TestPMCP
using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params, evaluate!
using CrystalShift: newton!, evaluate_residual!, BackgroundModel
using CrystalShift: PeakModCP, FixedPseudoVoigt, PhaseModel, full_optimize!
using LinearAlgebra
using CovarianceFunctions: EQ
using Random: rand
using Test
# using Plots

verbose = false
residual_tol = 0.1 # tolerance for residual norm after optimization
maxiter = 512 # appears to be required for phase combinations in particular

# Global
std_noise = 0.1
mean_θ = [1., .5, .2]
std_θ = [.05, 2., 1.]
# newton_lambda = 1e-2 TODO: make this passable to the newton optimization

test_path = "../data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]))
x = collect(8:.1:60)

pmcp = PeakModCP(CrystalPhase(String(s[1]), 0.1, FixedPseudoVoigt(0.01)),x, 10)
pmcp = PeakModCP(pmcp, (0.5.+0.5.*rand(10)).*get_free_params(pmcp))
y = zero(x)
evaluate!(y, pmcp, x)
# for i in eachindex(cs)
# c = optimize()
@time c = full_optimize!(PhaseModel(cs[1:1]), x, y, std_noise, mean_θ, std_θ;
                objective = "LS", method = LM,
                regularization = true,
                loop_num = 8,
                peak_shift_iter = 32,
                mod_peak_num = 32,
                peak_mod_mean = [1.],
                peak_mod_std = [.5],
                peak_mod_iter = 32,
                 verbose = false)
t = zero(x)
evaluate!(t, c, x)
@test norm(t-y) < 0.2

@time c = full_optimize!(cs[1], x, y, std_noise, mean_θ, std_θ;
                objective = "LS", method = LM,
                regularization = true,
                loop_num = 8,
                peak_shift_iter = 32,
                mod_peak_num = 32,
                peak_mod_mean = [1.],
                peak_mod_std = [.5],
                peak_mod_iter = 32, verbose = false)

t = zero(x)
evaluate!(t, c, x)
@test norm(t-y) < 0.2

@time c = full_optimize!(cs[1:1], x, y, std_noise, mean_θ, std_θ;
                objective = "LS", method = LM,
                regularization = true,
                loop_num = 8,
                peak_shift_iter = 32,
                mod_peak_num = 32,
                peak_mod_mean = [1.],
                peak_mod_std = [.5],
                peak_mod_iter = 32,
                 verbose = false)

t = zero(x)
evaluate!(t, c, x)
@test norm(t-y) < 0.2

evaluate_residual!(c, x, t)
@test norm(t) < 10^-10
# plot(x, y)
# plt = plot(x, y)
# plot!(x, evaluate!(zero(x), c, x))
# display(plt)

noise_intensity = 0.1
noise = noise_intensity.*(1 .+ sin.(0.2x))
noisy_data = noise .+ y
bg = BackgroundModel(x, EQ(), 10, rank_tol=1e-3)

@time c = full_optimize!(PhaseModel(cs[1:1], nothing, bg), x, noisy_data,
                std_noise, mean_θ, std_θ;
                objective = "LS", method = LM,
                regularization = true,
                loop_num = 8,
                peak_shift_iter = 32,
                mod_peak_num = 32,
                peak_mod_mean = [1.],
                peak_mod_std = [.5],
                peak_mod_iter = 32,
                verbose = false)

t = zero(x)
evaluate!(t, c, x)
@test norm(t-noisy_data) < 0.2

evaluate_residual!(c, x, t)
@test norm(t) < 10^-10
# plt = plot(x, noisy_data)
# plot!(x, evaluate!(zero(x), c, x))
# display(plt)
# end
# @time pm_cp = PeakModCP.(c.CPs, (x, ), (32, ))
# @time pm_cp = get_PeakModCP(c, x, 32)

# @time cp = optimize!(pm_cp, x, y, std_noise, peak_mean_θ, peak_std_θ;
#                   objective = objective, method = bfgs, maxiter = 16,
#                   regularization = true, verbose = false)
end