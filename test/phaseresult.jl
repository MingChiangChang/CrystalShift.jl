module Testphaseresult
using CrystalShift
using CrystalShift: PhaseResult, StripeResult, PseudoVoigt, get_center
using Test

verbose = false
residual_tol = 0.1 # tolerance for residual norm after optimization
maxiter = 1024 # appears to be required for phase combinations in particular

# Global
std_noise = 1e-3
mean_θ = [1., 1, .2]
std_θ = [.02, 1., 1.]
# newton_lambda = 1e-2 TODO: make this passable to the newton optimization

test_path = "../data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = @. CrystalPhase(String(s[1:end-1]), (0.1, ), (PseudoVoigt(0.5), ))
x = collect(8:.1:60)

a = PhaseResult([cs[1]], cs[1].name, 0.1, [1.,2.,3.], true)
b = PhaseResult([cs[2]], cs[2].name, 0.1, [1.,2.,3.], false)
c = PhaseResult([cs[3]], cs[3].name, 0.1, [1.,2.,3.], false)

sp = StripeResult([a,b,c], 0.5, 10, 10, 800, 800)
@test get_center(sp) == a
end