using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params
using CrystalShift: newton!, get_free_lattice_params
using CrystalShift: BackgroundModel, evaluate!, PhaseModel
using CovarianceFunctions
using CovarianceFunctions: EQ

using LinearAlgebra
using Random: rand
using Test
using Plots

verbose = false
residual_tol = 0.1 # tolerance for residual norm after optimization
maxiter = 512 # appears to be required for phase combinations in particular

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

noise_intensity = 0.1
noise = noise_intensity.*(1 .+ sin.(0.2x))
noisy_data = noise .+ data
bg = BackgroundModel(x, EQ(), 100)

PM = PhaseModel(cs[1:1], bg)
noisy_data = convert(Vector{Real}, noisy_data)

c = optimize!(PM, x, noisy_data, std_noise, mean_θ, std_θ,
            objective = "LS", method = LM, maxiter = maxiter,
            regularization = true, verbose = verbose)