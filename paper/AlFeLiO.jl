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
using JSON, Dates

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
# std_noise = 5e-2 #5e-3
# mean_θ = [1., .5, .2]
# std_θ = [0.05, 0.05, .05]
method = LM
opt_mode = Simple
objective = LeastSquares()
amorphous = false
K = 5
max_num_phases = 3
n_expand = 3
noise_level = 0.03
n_test = 231

std_noise = 1e-2
mean_θ = [1., .5, .5]
# std_θ = [0.05, .05, .05]
std_θ = [0.1, .05, .1]

###### Create CrystalPhase Objects #####
test_path = "data/AlFeLiO/sticks.csv"
open(test_path, "r") do f
    global cs = CrystalPhase(f, 0.1, Gauss() )
end

x = collect(15.0:.1:79.9)
d = npzread("data/AlFeLiO/alfeli.npy")
# d = npzread("alfeli_narrow.npy")
result_node = Vector{Vector{Node}}()

sol_path = "data/AlFeLiO/sol.csv"
sol = open(sol_path, "r")

t = split(read(sol, String), "\n")
gt = get_ground_truth(t)[1:n_test, :]
answer = zeros(Int64, (length(t), K, 7))

for i in tqdm(1:n_test)
    solution = split(t[i], ",")
    col = parse(Int, solution[1])
    # d[i, :] ./= maximum(d[i, :])
    # d[i, :] += abs.(noise_level .* randn(650)) # Add noise
    # d[i, :] ./= maximum(d[i, :])

    d[i, :] ./= maximum(d[i, :])
    # d[i, :] .*= 100
    # σₚ = sqrt.(((d[i, :])) .+ 1.)
    # σₚ = zero(sqrt.(d[i, :]))
    # d[i, :] += σₚ .* randn(650) # ADD POISSON NOISE
    # d[i, :] .= max.(d[i, :], 0)
    # σₚ ./= maximum(d[i, :])

    noise_level = [0.01]
    d[i, :] += abs.(noise_level[1] .* rand(650))

    d[i, :] ./= maximum(d[i, :])
    y = d[i, :]

    tree = Lazytree(cs, x)
    # Perform CrystalShift tree search
    result = search!(tree, x, y, max_num_phases, n_expand, amorphous, false, 5., std_noise, mean_θ, std_θ,
                        optimize_mode=opt_mode, method=method, maxiter=512, regularization=true, em_loop_num=1)
    if !amorphous
        result = result[2:end]
    end
    global result = vcat(result...)
    # Obtain the probabilities of the results
    global prob = get_probabilities(result, x, y, mean_θ, std_θ, objective=objective)

    highest = sortperm(prob, rev=true)[1:K]
    i_max = highest[1]
    plt = plot(x, y, label="Original", title="$(i)")
    plot!(x, result[i_max](x), label="Optimized")
    display(plt)
    push!(result_node, result[highest])
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

save_dict = Dict{String, Any}()
save_dict["std_noise"] = std_noise
save_dict["std_theta"] = std_θ
save_dict["mean_theta"] = mean_θ
save_dict["top_1_acc"] = top_k_accuracy(answer, gt, 1)
save_dict["top_5_acc"] = top_k_accuracy(answer, gt, 5)
save_dict["answer"] = answer
save_dict["gt"] = gt
save_dict["precision"] = precision(answer=answer[1:n_test,1,:], ground_truth=gt, verbose=false)
save_dict["recll"] = recall(answer=answer[1:n_test,1,:], ground_truth=gt, verbose=false)

open("alfelio_$(Dates.format(now(), "yyyy-mm-dd_HH:MM")).json", "w") do f
    JSON.print(f, save_dict)
end