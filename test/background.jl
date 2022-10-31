module TestBackground
using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params
using CrystalShift: newton!, get_free_lattice_params
using CrystalShift: BackgroundModel, evaluate!, PhaseModel, Lorentz, PseudoVoigt
using CrystalShift: _prior, get_param_nums, lm_prior!
using CovarianceFunctions
using CovarianceFunctions: EQ

using LinearAlgebra
using Random: rand
using Test

verbose = false
residual_tol = 0.1 # tolerance for residual norm after optimization
maxiter = 2000 # appears to be required for phase combinations in particular

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


function synthesize_data(cp::CrystalPhase, x::AbstractVector)
    params = get_free_lattice_params(cp)
    interval_size = 0.025
    scaling = (interval_size.*rand(size(params, 1),) .- interval_size/2) .+ 1
    @. params = params*scaling
    params = [params..., 1., 0.2, 0.5]
    r = evaluate(cp, params, x)
    normalization = maximum(r)
    params[end-1] /= normalization
    if verbose
        println("synthesize_data: ", params)
    end

    return r/normalization, params
end

data, params = synthesize_data(cs[1], x)

# TODO: Test pure background fit
noise_intensity = 0.1
noise = noise_intensity.*(1 .+ sin.(0.2x))
noisy_data = noise .+ data
bg = BackgroundModel(x, EQ(), 10, rank_tol=1e-4)

PM = PhaseModel(cs[1:1], nothing, bg)
noisy_data = convert(Vector{Real}, noisy_data)

@time c = optimize!(PM, x, noisy_data, std_noise, mean_θ, std_θ,
            objective = "LS", method = LM, optimize_mode = Simple, maxiter = maxiter,
            regularization = true, verbose = verbose)

y = zero(x)

# using Plots
# plt = plot(x, noisy_data, label="Ground truth")
# plot!(x, evaluate!(zero(x), c, x), label="Fitted graph")
# display(plt)
@test norm(evaluate!(y, c, x) .- noisy_data) < 0.1

evaluate(PM.background, x)
evaluate(PM.background, rand(get_param_nums(PM.background)), x)
_prior(PM.background, rand(get_param_nums(PM.background)))
lm_prior!(rand(get_param_nums(PM.background)), PM.background)

end