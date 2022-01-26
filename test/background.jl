using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params
using CrystalShift: newton!, get_free_lattice_params

using LinearAlgebra
using Random: rand
using Test

verbose = false
residual_tol = 0.1 # tolerance for residual norm after optimization
maxiter = 128 # appears to be required for phase combinations in particular

# Global
std_noise = 1e-3
mean_θ = [1., 1, .2]
std_θ = [.02, 1., 1.]
# newton_lambda = 1e-2 TODO: make this passable to the newton optimization

test_path = "data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]))
x = collect(8:.1:60)