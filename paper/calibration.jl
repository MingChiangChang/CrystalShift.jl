using CrystalShift
using CrystalTree
using CrystalTree: Lazytree, search!, approximate_negative_log_evidence, get_phase_ids, LeastSquares, get_probabilities
using CrystalShift: Lorentz, get_free_lattice_params, extend_priors, get_free_params
using Combinatorics
using ProgressBars
using Measurements
using NPZ

using StatsBase
using Plots
using Statistics: cor

"""
    get_bin(prob::Any)

    Get the bin number for a given probability value.
"""
function get_bin(prob)
    if isnan(prob)
        p = 1
    else
        p = Int64(floor(prob*10))+1
    end
    if p > 10
        p = 10
    end
    # println("$(prob) $(p)")
    return p
end

function get_mod_phase_ids(pm)
    ids = get_phase_ids(pm)
    for i in eachindex(ids)
        ids[i] += 1
    end
    Set(ids)
end

std_noises = .05
mean_θ = [1., .5, .1]
std_θ = [.05, .05, .05]
amorphous = false
k = 2
runs = 1000
correct_count = 0

test_path = "paper/data/calibration/sticks.csv"
data_path = "paper/data/calibration/calibration_data_normal_0.03_easy.npy"
test_data = npzread(data_path)


x = collect(8:.1:40) # Q values
totl = zeros(Int64, 10)
correct = zero(totl)
totl_prob = zeros(Float64, 10)
phase_correct = zeros(Int64, 10)
phase_totl = zeros(Int64, 10)

for i in tqdm(1:runs)
    #### Create CrystalPhase Object ####
    open(test_path, "r") do f
        global cs = CrystalPhase(f, 0.1, Lorentz())
    end
    y = test_data[i, 1:end-2] # The last two columns are the ground truth phase ids
    test_comb = test_data[i, end-1:end]
    if test_comb[end] == 0.0
        pop!(test_comb)
    end
    test_comb = Set(test_comb)

    LT = Lazytree(cs, x)
    # Where the CrystalShift perform search
    results = search!(LT, x, y, 2, k, amorphous, false, 5., std_noise, mean_θ, std_θ,
                    method=LM, objective="LS", optimize_mode=EM,
                    maxiter=256, em_loop_num=3,
                    regularization=true)
    if !amorphous
        results = results[2:end]
    end
    results = reduce(vcat, results)
    # Obtain the probabilities of the results
    probs = get_probabilities(results, x, y, mean_θ, std_θ, renormalize=true, normalization_constant=.5)

    prob_of_phase = zeros(Float64, 5)
    for j in eachindex(results)
        for k in eachindex(results[j].phase_model.CPs)
            ind = results[j].phase_model.CPs[k].id + 1
            prob_of_phase[ind] += probs[j]
        end
    end

    for j in eachindex(prob_of_phase)
        ind = get_bin(prob_of_phase[j])
        phase_totl[ind] += 1
        if j in test_comb
            phase_correct[ind] += 1
        end
    end

    ind = argmax(probs)
    ss = Set([results[ind].phase_model.CPs[i].id+1 for i in eachindex(results[ind].phase_model.CPs)])
    answer = test_comb
    ss == answer && (global correct_count += 1)

    for j in eachindex(results)
        bin_num = get_bin(probs[j])
        totl[bin_num] += 1
        totl_prob[bin_num] += probs[j]
        if get_mod_phase_ids(results[j]) == test_comb
            correct[bin_num] += 1
        end
    end
end


pearson = cor([0.05+0.1*i for i in 0:9], correct./totl)

plt = plot([0., 1.], [0., 1.],
        linestyle=:dash, color=:black,
        legend=false, figsize=(10,10), dpi=300,
        xlims=(0, 1), ylims=(0, 1), xtickfontsize=10, ytickfontsize=10,
        xlabelfontsize=12, ylabelfontsize=12, markersize=5,
        title="k=$(k)\nstd_noise=$(std_noise), mean=$(mean_θ)\n std=$(std_θ)\n runs=$(runs) pearson=$(pearson)\n accuracy=$(correct_count/runs)")

calibration = correct ./ totl

for i in eachindex(calibration)
    if isnan(calibration[i])
        calibration[i] = 0
    end
end

plot!(collect(0.05:.1: 0.95), calibration)
scatter!(collect(0.05:.1: 0.95), calibration)
plot!(totl_prob ./ totl, calibration)
scatter!(totl_prob ./ totl, calibration)

font(20)
xlabel!("Predicted probabilities")
ylabel!("Frequency of correct matches")
display(plt)

pearson = cor([0.05+0.1*i for i in 0:9], phase_correct./phase_totl)

plt = plot([0., 1.], [0., 1.],
        linestyle=:dash, color=:black,
        legend=false, figsize=(10,10), dpi=300,
        xlims=(0, 1), ylims=(0, 1), xtickfontsize=10, ytickfontsize=10,
        xlabelfontsize=12, ylabelfontsize=12, markersize=5,
        title="k=$(k)\nstd_noise=$(std_noise), mean=$(mean_θ)\n std=$(std_θ)\n runs=$(runs) pearson=$(pearson)\n accuracy=$(correct_count/runs)")

calibration = phase_correct ./ phase_totl

for i in eachindex(calibration)
    if isnan(calibration[i])
        calibration[i] = 0
    end
end

plot!(collect(0.05:.1: 0.95), calibration)
scatter!(collect(0.05:.1: 0.95), calibration)

plot!(totl_prob ./ totl, calibration)
scatter!(totl_prob ./ totl, calibration)

font(20)
xlabel!("Predicted phase probabilities")
ylabel!("Frequency of correct phase matches")
display(plt)

t = Dict{Any, Any}()
t["std_noise"] = std_noise
t["mean_theta"] = mean_θ
t["std_theta"] = std_θ
t["runs"] = runs
t["accuracy"] = correct_count/runs
t["totl"] = totl
t["correct"] = correct
t["phase_correct"] = phase_correct
t["totl_prob"] = totl_prob
t["phase_totl"] = phase_totl

using JSON
using Dates

open("BTO_Gaussian_=0.05_easy_test_$(Dates.format(now(), "yyyy-mm-dd_HH:MM")).json", "w") do f
    JSON.print(f, t)
end