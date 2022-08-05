module testfixedbackground
using CrystalShift
using CrystalShift: CrystalPhase, optimize!, evaluate, get_free_params
using CrystalShift: newton!, get_free_lattice_params
using CrystalShift: BackgroundModel, evaluate!, PhaseModel, Lorentz, PseudoVoigt
using CrystalShift: _prior, get_param_nums, lm_prior!
using CrystalShift: FixedBackground

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

# construct polynomial background

poly_bg = zero(x)
@. poly_bg = (x-2) * (x-20) * (x-30) * (x-50) * (x-70)
poly_bg .-= minimum(poly_bg)
poly_bg ./= maximum(poly_bg) * 2

data += poly_bg
# Add noise
noise_intensity = 0.05
noise = noise_intensity.*rand(size(x, 1))
noisy_data = noise .+ data

noisy_data ./= maximum(noisy_data)

basis = zero(x)
@. basis = (x-2) * (x-20) * (x-30) * (x-50) * (x-70)
basis .-= minimum(basis)
basis ./= maximum(basis) * 2

FB = FixedBackground(basis, 1., 0., 0.)
pm = PhaseModel([cs[1]], nothing, FB)
c = optimize!(pm, x, noisy_data, std_noise, mean_θ, std_θ,
            objective = "LS", method = LM, maxiter = maxiter,
            regularization = true, verbose = verbose)

# Check MSE
@test (norm(noisy_data .- evaluate!(zero(x), c, x))^2 / 521) < 0.01

end