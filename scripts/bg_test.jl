using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params
using CrystalShift: newton!, get_free_lattice_params, PseudoVoigt
using CrystalShift: BackgroundModel, evaluate!, PhaseModel, Lorentz, FixedPseudoVoigt
using CovarianceFunctions
using CovarianceFunctions: EQ

using LinearAlgebra
using Random: rand
using Test
using HDF5
using Plots

verbose = false
residual_tol = 0.1 # tolerance for residual norm after optimization
maxiter = 512 # appears to be required for phase combinations in particular

# Global
std_noise = .1
mean_θ = [1., 1., .2]
std_θ = [1., .5, .5]
# newton_lambda = 1e-2 TODO: make this passable to the newton optimization

test_path = "/Users/ming/Downloads/sticks_1.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = @. CrystalPhase(String(s[1:end-1]), (.1, ), (PseudoVoigt(0.8), ))

h5 = h5open("/Users/ming/Downloads/31_21F66_SnOx_all_1d.h5")
t = h5["exp"]["tau_5866_T_1200"]["125"]["integrated_1d"]
x = t[100:900,1]
y = t[100:900,2]
y ./= maximum(y)

bg = BackgroundModel(x, EQ(), 10)
PM = PhaseModel([cs[3], cs[7]], bg)

# plot(x, y)

c = optimize!(PM, x, y, std_noise, mean_θ, std_θ,
            objective = "LS", method = LM, maxiter = maxiter,
            regularization = true, verbose = verbose)

plt = plot(x, y)
plot!(x, evaluate!(zero(x), c, x))
display(plt)