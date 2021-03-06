using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params, Lorentz
using CrystalShift: newton!, get_free_lattice_params, get_fraction, PseudoVoigt
using CrystalShift: evaluate!,get_param_nums, FixedPseudoVoigt

using LinearAlgebra
using Random: rand
using BenchmarkTools
using Plots

verbose = false
residual_tol = 0.1 # tolerance for residual norm after optimization
maxiter = 512 # appears to be required for phase combinations in particular

# Global
std_noise = 1.
mean_θ = [1., .5, .2]
std_θ = [.05, 2., 1.]
# newton_lambda = 1e-2 TODO: make this passable to the newton optimization

function synthesize_data(cp::CrystalPhase, x::AbstractVector)
    params = get_free_lattice_params(cp)
    interval_size = 0.025
    scaling = (interval_size.*rand(size(params, 1),) .- interval_size/2) .+ 1
    @. params = params*scaling
    params = [params..., 1., 0.2, 0.5*rand(1)...]
    r = zero(x)
    evaluate!(r, cp, params, x)
    normalization = maximum(r)
    params[end-2] /= normalization
    return r/normalization, params
end

test_path = "data/Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]), (0.1,), (Lorentz(),))
x = collect(8:.1:60)
y = zero(x)

test, params = synthesize_data(cs[1], x)
# test = evaluate!(zero(x), cs[1], params, x)
# @benchmark c = optimize!(cs[1], x, test, std_noise, mean_θ, std_θ;
#                   method =LM,
#                   maxiter = maxiter, regularization = true, verbose=false)

@benchmark evaluate!(zero(x), cs[1], x)
# plt = plot(x, test)
# plot!(x, evaluate!(y, c, x))
# display(plt)

# t = rand(get_param_nums(cs[1]))
# tt = rand(get_param_nums(cs[1:3]))
# println("single phase")
# @btime evaluate!(y, cs[1], t, x)
# println("3 phases")
# @btime evaluate!(y, cs[1:3], tt, x)
# println("")