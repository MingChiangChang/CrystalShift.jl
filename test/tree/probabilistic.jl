module TestProbabilistic
using Test
using CrystalShift: approximate_negative_log_evidence, get_probabilities
using CrystalShift: bestfirstsearch, find_first_unassigned, get_phase_ids, res_bfs
using CrystalShift
using CrystalShift: CrystalPhase, optimize!, get_free_params, Lorentz
using CrystalShift: extend_priors, evaluate_residual!, evaluate!
using ProgressBars
using LinearAlgebra


std_noise = 1e-2
mean_θ = [1., 1., .1]
std_θ = [0.05, 5., .1]

# CrystalPhas object creation
path = "../data/"
phase_path = path * "tree_test_sticks.csv"
f = open(phase_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

if s[end] == ""
    pop!(s)
end

cs = Vector{CrystalPhase}(undef, size(s))
@. cs = CrystalPhase(String(s), (0.1,), (Lorentz(),))
max_num_phases = 3
tree = Tree(cs[1:15], max_num_phases)
x = collect(8:.035:45)
y = zero(x)

for node in tree.nodes[2:3]
    node.phase_model(x, y)
end

noise = rand(size(x, 1))


y /= maximum(y)  # normalize targets

method = LM
objective = CrystalShift.LeastSquares()
results = res_bfs(tree, x, y, std_noise, mean_θ, std_θ, 15,
                        method=method, objective=objective, optimize_mode=Simple,
                        maxiter=1000, regularization=true)

num_nodes = find_first_unassigned(results) - 1
results = results[1:num_nodes]
neg_log_prob = zeros(num_nodes)
for (i, result) in enumerate(results)
    θ = get_free_params(result.phase_model)
    # orig = [p.origin_cl for p in result[i].phase_model]
    full_mean_θ, full_std_θ = extend_priors(mean_θ, std_θ, result.phase_model.CPs)
    neg_log_prob[i] = approximate_negative_log_evidence(result, θ, x, y,
                        std_noise, full_mean_θ, full_std_θ, objective)
end

i_min = argmin(neg_log_prob)

true_phase_ids = Set((0, 1))
phase_ids = Set(get_phase_ids(results[i_min]))
@test phase_ids == true_phase_ids

# test for get_probabilities
prob = get_probabilities(results, x, y, std_noise, mean_θ, std_θ)

# validity as probability distribution
@test sum(prob) ≈ 1
@test all(≥(0), prob)

# non-zero probabilities are in reverse sort order of negative log probabilities
ind = @. prob != 0
@test sortperm(prob[ind]) == reverse(sortperm(neg_log_prob[ind]))

# highest-probability phase is correct
i_min = argmax(prob)
true_phase_ids = Set((0, 1))
phase_ids = Set(get_phase_ids(results[i_min]))
@test phase_ids == true_phase_ids

# testing without renormalization
prob = get_probabilities(results, x, y, std_noise, mean_θ, std_θ, renormalize=false)
@test sum(prob) ≈ 1
@test all(≥(0), prob)
@test exp.(-neg_log_prob) / sum(exp.(-neg_log_prob)) ≈ prob

# Hessians for get_probabilities: exact (block JᵀJ + hand-written curvature) and Gauss–Newton
# (2JᵀJ + prior), against ForwardDiff
using ForwardDiff
using CrystalShift: least_squares_hessian, ls_regularizer, ls_objective, _weighted_residual!, BackgroundModel,
                    EQ, PhaseModel, get_param_nums, Gauss, PseudoVoigt, FixedPseudoVoigt, USE_ANALYTIC_JACOBIAN
function log_point(pm)
    θ = get_free_params(pm)
    m, s = extend_priors(mean_θ, std_θ, pm.CPs)
    nc = get_param_nums(pm.CPs)
    log_θ = copy(θ); log_θ[1:nc] .= log.(log_θ[1:nc])
    log_θ, t -> ls_regularizer(pm, t, log.(m), s), log.(m), s
end
weighted_r(pm, t) = _weighted_residual!(pm, t, x, y, zero(y), zeros(eltype(t), length(y)), std_noise)

# fitted nodes (Lorentz), every supported profile away from the optimum, and a background model
hessian_models = [[r.phase_model for r in results[1:min(end, 5)]];
                  [PhaseModel(CrystalPhase.(String.(s[[1, 6]]), (0.1,), (P,))) for P in (Gauss(), PseudoVoigt(0.5), FixedPseudoVoigt(0.5))];
                  PhaseModel(cs[[1]], nothing, BackgroundModel(x, EQ(), 8., 10., rank_tol = 1e-3))]

@testset "Exact Hessian matches ForwardDiff" begin
    for pm in hessian_models, analytic in (true, false)
        log_θ, R, mlog, sd = log_point(pm)
        USE_ANALYTIC_JACOBIAN[] = analytic
        H = least_squares_hessian(pm, log_θ, x, y, zero(y), std_noise, R; exact = true)
        USE_ANALYTIC_JACOBIAN[] = true
        @test H ≈ ForwardDiff.hessian(t -> ls_objective(pm, t, x, y, zero(y), std_noise, mlog, sd), log_θ) rtol = 1e-10
    end
end

@testset "Gauss–Newton Hessian matches 2JᵀJ + prior Hessian" begin
    for pm in hessian_models
        log_θ, R, _, _ = log_point(pm)
        J = ForwardDiff.jacobian(t -> weighted_r(pm, t), log_θ)
        H = least_squares_hessian(pm, log_θ, x, y, zero(y), std_noise, R; exact = false)
        @test H ≈ 2J'J + ForwardDiff.hessian(R, log_θ) rtol = 1e-10
    end
end

@testset "get_probabilities hessian options" begin
    prob_fd = get_probabilities(results, x, y, std_noise, mean_θ, std_θ, hessian = :forwarddiff)
    @test get_probabilities(results, x, y, std_noise, mean_θ, std_θ) ≈ prob_fd # :exact is the default
    prob_gn = get_probabilities(results, x, y, std_noise, mean_θ, std_θ, hessian = :gauss_newton)
    @test sum(prob_gn) ≈ 1
    @test all(≥(0), prob_gn)
    @test argmax(prob_gn) == argmax(prob_fd)
    θ = get_free_params(results[1].phase_model)
    m, sd = extend_priors(mean_θ, std_θ, results[1].phase_model.CPs)
    @test_throws ArgumentError approximate_negative_log_evidence(results[1], θ, x, y, zero(y), std_noise, m, sd,
                                                                 objective; hessian = :newton)
end
println("End of probabilistic.jl test")
end # TestProbabilistic module
