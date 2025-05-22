# Ta-Sn-O Phase diagram and lattice parameter
using NPZ
using Plots
using YAML: load_file

# using PhaseMapping: xray
using BackgroundSubtraction: mcbl
using LinearAlgebra
using JSON
using ProgressBars
using CSV
using DataFrames
using CrystalShift
using CrystalTree
using CrystalShift: parse_cond, PhaseResult, get_weighted_center, StripeResult, parse_cond, CrystalPhase
using CrystalShift: Lorentz, get_free_params, extend_priors
using CrystalTree: Lazytree, search!, approximate_negative_log_evidence, LeastSquares, get_probabilities

# Overloading
using JSON: lower
include("nmf.jl")
JSON.lower(::Lorentz) = nothing

sol = load_file("data/TaSnO/gt.yaml")

########## Constant Declaration ##########
const RANK = 4
const THRESH = 0.3

std_noise = .03
mean_θ = [1., .5, .2]
std_θ = [0.05,0.05,.05]
maxiter = 256
max_depth = 3
k = 3
amorphous = false
background = false
background_length = 5.

########## Load data ##########
dl = "data/TaSnO/"
# data = npzread(dl * "TaSnO_data.npy")
Ws = npzread(dl * "TaSnO_Ws.npy")
Hs = npzread(dl * "TaSnO_Hs.npy")
Ks = npzread(dl * "TaSnO_k.npy")
q = npzread(dl * "TaSnO_Q.npy")
uncer = npzread(dl * "TaSnO_uncer.npy")

open(dl * "TaSnO_conds.json", "r") do f
    global conditions = JSON.parse(f)
end

conds = parse_cond.(conditions, Float64)

df = DataFrame(CSV.File(dl * "TaSnO_DwellTpeak.csv"))
ta = df[!, "Ta.nmol_offgrid"]
sn = df[!, "Sn.nmol_offgrid"]
cation = zero(ta)
@. cation = ta/(ta+sn)

########## CrystalPhas object creation ##########
phase_path = dl * "TaSnO_sticks.csv"
open(phase_path, "r") do f
    global cs = CrystalPhase(f, 0.1, Lorentz())
end

########## Result Container Init ##########
wafer_result = Vector{StripeResult}(undef, size(Ws)[1])
totl = 0 # total ground truth phases
correct = 0 # Absolute correct cases
tp = 0 # True positive
fp = 0 # False positive
fn = 0 # False negative

########## CrystalShift doing work ##########
for i in tqdm(1:size(Ws, 1))
    condition = conds[i]

    BW = Ws[i, :, :]'
    BH = Hs[i, :, :]
    BK = Ks[i, :]

    stripe = Vector{Vector{PhaseResult}}(undef, 4)
    center = convert(Int64, round(get_weighted_center(BH))) # Weighted center with the activation

    isCenter = BitArray(undef, RANK) # Indicator of whether this basis exist at the center
    for k in 1:RANK
        isCenter[k] = BH[k, center] > THRESH
    end

    for j in 1:RANK
        gt_names = sol[i][j]
        nc = maximum(BW[:, j])
        BW[:, j] ./= nc # maximum(BW[:, j])
        b = mcbl(BW[:, j], q[i,:], 7) # Background prediction
        y_bg_sub = BW[:, j] - b
        y_uncer = uncer[i, BK[j], :] ./ nc
        if "a" in gt_names
            stripe[j] =  [ PhaseResult(cs[1], "non-crystalline", BH[j, :], y_bg_sub, isCenter[j])] # use cs[1] as dummy
            continue
        end
        @. y_bg_sub = max(y_bg_sub, 0)

        tree = Lazytree(cs, q[i,:])
        # Perform tree search
        results = search!(tree, q[i,:], y_bg_sub, y_uncer, max_depth, k, amorphous, background, background_length,
                            std_noise, mean_θ, std_θ,
                            method=LM, objective="LS", maxiter=maxiter, optimize_mode=EM, em_loop_num=5,
                            regularization=true)

        results = results[2:end]
        results = reduce(vcat, results)

        # obtain probabilistic labeling results
        probs = get_probabilities(results, q[i,:], y_bg_sub, mean_θ, std_θ)

        result_node = results[argmax(probs)]
        ########## Recording Metrics ###########
        gt_names = sol[i][j]
        phase_names = getproperty.(result_node.phase_model.CPs, :name)
        phase_names = map(x->split(x, '_')[1], phase_names)

        phase_names = String.(phase_names)
        global totl += 1
        if phase_names == sol[i][j]
            global correct += 1
        end

        for phase_name in phase_names
            if phase_name in sol[i][j]
                global tp += 1
            else
                global fp += 1
            end
        end

        for gt_phase_name in sol[i][j]
            if gt_phase_name ∉ phase_names
                global fn += 1
            end
        end

        ########## Plotting ##########
        plt = plot(q[i,:], y_bg_sub, xtickfontsize=10, ytickfontsize=10, lw=4, label="Diffraction")
        n = ""
        # for p in result_node.phase_model
        #     n = n * p.name * "\n"
        # end
        for p in result_node.phase_model
            plot!(q[i, :], p.(q[i, :]), label=p.name, ylims=(0.0, 1.0), lw=2)
        end
        title!("$(condition) $(j)")
        xlabel!("Q (1/nm)", fontsiz=10)
        ylabel!("a.u.", fintsize=10)
        # savefig("figures/$(condition) $(j).png")
        display(plt)

        stripe[j] =  [ PhaseResult(p, p.name, BH[j, :], y_bg_sub, isCenter[j])
                    for p in result_node.phase_model.CPs]
    end
    subdf = df[(df.Tpeak .== conds[i][1]) .& (df.dwell .== conds[i][2]), :]
    local ta = subdf[!, "Ta.nmol_offgrid"]
    local sn = subdf[!, "Sn.nmol_offgrid"]
    local cation = ta./(ta.+sn)
    wafer_result[i] =  StripeResult(stripe, cation, subdf.xcenter[1], subdf.ycenter[1], conds[i]...)
end

# Store data
using Dates
open("data/TaSnO_metric_$(Dates.format(now(), "yyyy-mm-dd_HH:MM")).json", "w") do f
    JSON.print(f, wafer_result)
end

metrics = Dict{String, Any}()
metrics["std_noise"] = std_noise
metrics["mean_theta"] = mean_θ
metrics["std_theta"] = std_θ
metrics["max_depth"] = max_depth
metrics["k"] = k
metrics["maxiter"] = maxiter
metrics["totl"] = totl
metrics["correct"] = correct
metrics["tp"] = tp
metrics["fp"] = fp
metrics["fn"] = fn


open("TaSnO_$(Dates.format(now(), "yyyy-mm-dd_HH:MM")).json", "w") do f
    JSON.print(f, metrics)
end