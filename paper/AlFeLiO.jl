using ProgressBars

using CrystalTree
using CrystalTree: approximate_negative_log_evidence, find_first_unassigned
using CrystalTree: Lazytree, search_k2n!, search!, cast, LeastSquares
using CrystalTree: get_phase_number, get_ground_truth, precision, recall
using CrystalTree: in_top_k, top_k_accuracy, get_probabilities
using CrystalShift
using CrystalShift: get_free_params, extend_priors, Lorentz, evaluate_residual!, PseudoVoigt
using CrystalShift: Gauss, FixedPseudoVoigt

using Plots
using LinearAlgebra
using Base.Threads
using DelimitedFiles
using NPZ

function node_under_improvement_constraint(nodes, improvement, x, y)
    min_res = (Inf, 1)
    for i in eachindex(nodes)
        copy_y = copy(y)
        res = norm(evaluate_residual!(nodes[i].phase_model, x, copy_y))
        if min_res[1] - res > improvement
            min_res = (res, i)
        end
    end
    return  nodes[min_res[2]]
end

####### CrystalShift setup #######
std_noise = 5e-3
mean_θ = [1., .5, .5]
std_θ = [0.05, 1.0, .05]
method = LM
opt_mode = Simple
objective = LeastSquares()
amorphous = false
K = 5
max_num_phases = 3
noise_level = 0.03


test_path = "paper/data/AlFeLiO/sticks.csv"
open(test_path, "r") do f
    global cs = CrystalPhase(f, 0.1, Gauss() )
end

x = collect(15.0:.1:79.9)
d = npzread("paper/data/AlFeLiO/alfeli.npy")
result_node = Vector{Vector{Node}}()

sol_path = "paper/data/AlFeLiO/sol.txt"
sol = open(sol_path, "r")

t = split(read(sol, String), "\n")
gt = get_ground_truth(t)
answer = zeros(Int64, (length(t), K, 7))

for i in tqdm(eachindex(t))
    solution = split(t[i], ",")
    col = parse(Int, solution[1])
    d[i, :] ./= maximum(d[i, :])
    d[i, :] += abs.(noise_level .* randn(650)) # Add noise
    d[i, :] ./= maximum(d[i, :])
    y = d[i, :]

    tree = Lazytree(cs, x)
    result = search!(tree, x, y, 3, 3, 1., amorphous, false, 5., std_noise, mean_θ, std_θ,
                        maxiter=512, regularization=true)
    if !amorphous
        result = result[2:end]
    end
    global result = vcat(result...)

    global prob = Vector{Float64}(undef, length(result))
    @threads for i in eachindex(result)
        θ = get_free_params(result[i].phase_model)
        full_mean_θ, full_std_θ = extend_priors(mean_θ, std_θ, result[i].phase_model.CPs)
        prob[i] = approximate_negative_log_evidence(result[i], θ, x, y, std_noise, full_mean_θ, full_std_θ, objective, true) # put std_noise in for fix σₙ
    end

    lowest = sortperm(prob)[1:K]
    i_min = lowest[1]
    plt = plot(x, y, label="Original", title="$(i)")
    plot!(x, result[i_min](x), label="Optimized")
    display(plt)
    push!(result_node, result[lowest])
end

for i in eachindex(result_node)
    one_phase_answer = zeros(Int64, (K, 7))
    for j in eachindex(result_node[i])
        re = zeros(Int64, 7)
        for k in eachindex(result_node[i][j].phase_model.CPs)
            re[get_phase_number(result_node[i][j].phase_model.CPs[k].name)] += 1
        end
        one_phase_answer[j, :] = re
    end
    answer[i, :, :] = one_phase_answer
end

using JSON, Dates
d = Dict{String, Any}()
d["std_noise"] = std_noise
d["std_theta"] = std_θ
d["mean_theta"] = mean_θ
d["top_1_acc"] = top_k_accuracy(answer, gt, 1)
d["top_5_acc"] = top_k_accuracy(answer, gt, 5)
d["answer"] = answer
d["gt"] = gt
d["precision"] = precision(answer=answer[:,1,:], ground_truth=gt, verbose=false)
d["recll"] = recall(answer=answer[:,1,:], ground_truth=gt, verbose=false)

open("alfelio_$(Dates.format(now(), "yyyy-mm-dd_HH:MM")).json", "w") do f
    JSON.print(f, d)
end