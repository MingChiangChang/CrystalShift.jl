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
println("End of probabilistic.jl test")
end # TestProbabilistic module
