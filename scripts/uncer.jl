using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params, evaluate!, Lorentz
using CrystalShift: newton!, get_free_lattice_params, get_fraction, optimize_with_uncertainty!

using CovarianceFunctions: EQ
using LinearAlgebra
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

test_path = "data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]), (0.1,), (Lorentz(),))
x = collect(8:.1:60)
bg = BackgroundModel(x, EQ(), 10., 100., rank_tol=1.)

y = zero(x)
noise = 0.1rand(length(x))
y += noise
evaluate!(y, CrystalPhase(cs[1], [17.0, 4.86, 5.55, 1.58, 1.0, 0.2]), x)
pm, uncer = optimize_with_uncertainty!(PhaseModel(cs[1:1], nothing, bg), x, y, std_noise, mean_θ, std_θ;
                            method = LM,
                            maxiter = maxiter,
                            regularization = true,
                            verbose=false)
plot(x, y)