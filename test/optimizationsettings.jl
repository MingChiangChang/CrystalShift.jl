module Testoptsetting
using CrystalShift
using CrystalShift: extend_priors, PseudoVoigt, PeakModCP, OptimizationSettings
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


@test extend_priors(mean_θ, std_θ, cs[1:2]) == ([17.11299, 4.872, 5.548, 1.5810937693816631, 1.0, 0.2, 0.5, 4.072593, 1.0, 0.2, 0.5], [0.3422598, 0.09744, 0.11096, 0.031621875387633266, 1.0, 1.0, 10.0, 0.08145186000000001, 1.0, 1.0, 10.0])
pmcp = PeakModCP(cs[1], x, 10)
pmcp2 = PeakModCP(cs[2], x, 10)
@test extend_priors([1.], [.02], [pmcp, pmcp2]) == ([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02])

OptimizationSettings{Float64}(cs[1:2], std_noise, mean_θ,
           std_θ, 128, true, LM, "LS", Simple, false, 1e-6)
end