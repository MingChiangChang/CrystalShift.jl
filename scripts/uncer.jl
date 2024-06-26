using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params, evaluate!, Lorentz
using CrystalShift: newton!, get_free_lattice_params, get_fraction, optimize_with_uncertainty!
using CrystalShift: var_lognormal, get_eight_params, uncertainty, get_fraction_with_uncertainty
using CovarianceFunctions: EQ
using LinearAlgebra
using Random: rand
using Test
using Plots

using Measurements
# using Plots

verbose = false
residual_tol = 0.1 # tolerance for residual norm after optimization
maxiter = 2000 # appears to be required for phase combinations in particular

# Global
std_noise = .1
mean_θ = [1., .5, .2]
std_θ = [.05, 2., 1.]
# newton_lambda = 1e-2 TODO: make this passable to the newton optimization

test_path = "data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]), (0.1,), (Lorentz(),))
x = collect(8:.1:60)
# bg = BackgroundModel(x, EQ(), 10., 100., rank_tol=1.)

y = zero(x)
noise = 0.2rand(length(x))
y += noise
# evaluate!(y, CrystalPhase(cs[1], [17.0, 4.86, 5.55, 1.58, 1.0, 0.2]), x)
evaluate!(y, cs[1], x)
evaluate!(y, cs[2], x)
pm, uncer = optimize!(PhaseModel(cs[1:2], nothing, nothing), x, y, std_noise, mean_θ, std_θ;
                            method = LM, optimize_mode = WithUncer,
                            maxiter = maxiter,
                            regularization = true,
                            verbose = false)
plt = plot(x, y)
plot!(x, evaluate!(zero(x), pm, x))
display(plt)
println(uncer)
opt_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ)
#println(get_fraction_with_uncertainty(pm.CPs, x, y, false))
uncer = uncertainty(pm.CPs, x, y, opt_stn)
mean = log.(reduce(vcat, get_eight_params.(pm.CPs)))
std = sqrt.(var_lognormal.(mean, sqrt.(uncer)))
v = mean .± std

# exp.(mean .± std)
exp.(v)