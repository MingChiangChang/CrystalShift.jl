# Figures for the NMF + probabilistic phase determination work flow

using CrystalShift
using CrystalShift: Lorentz, get_free_params, extend_priors, FixedPseudoVoigt, Background
using CrystalShift: Crystal, OptionalPhases, get_param_nums, get_fraction
using CrystalTree
using CrystalTree: Lazytree, approximate_negative_log_evidence, LeastSquares
using BackgroundSubtraction
using CovarianceFunctions: EQ

using NPZ
using Plots
using LinearAlgebra
using ProgressBars
using Base.Threads

include("nmf.jl")

function get_phase_model_with_phase_names(phase_names::AbstractSet,
                                        phases::AbstractArray,
                                        background::Background=nothing,
                                        wildcard::OptionalPhases=nothing)
    l = []
    if isempty(phase_names)
        return PhaseModel(phases[l], wildcard, background)
    end
    for phase_name in phase_names
        for i in eachindex(phases)
            if phases[i].name == phase_name
                push!(l, i)
            end
        end
    end
    PhaseModel(phases[l], wildcard, background)
end

param_name = Dict{Int64, AbstractArray}()
param_name[1] = ["a"]
param_name[2] = ["a", "c"]
param_name[3] = ["a", "b", "c"]
param_name[4] = ["a", "b", "c", "β"]
param_name[6] = ["a", "b", "c", "α", "β", "γ"]

data = npzread("paper/data/TaSnO/TaSnO_data.npy")
q = npzread("paper/data/TaSnO/TaSnO_Q.npy")

##### Parameters #####
std_noise = .05
mean_θ = [1., .5, .1]
std_θ = [.001, 2., 1.]
RANK = 4
STRIPE_IDX = 336

##### Plotting heatmap #####
# heatmap(LinRange(-.75, .75, 76), q[336,50:end-50], transpose(data[336, 77:152,50:end-50]),
#         yflip=true, clim=(0, 10), xlabel="Relative Location (mm)", ylabel="q (nm⁻¹)",
#         ylabelfontsize=12, xlabelfontsize=12, xtickfontsize=10, ytickfontsize=10, colorbar_title="Intensity (a.u.)")


q = q[1, 1:900]
W, H, K = xray(Array(transpose(data[STRIPE_IDX, :, 1:900])), RANK)

##### Plotting bases #####
# plt = plot(xlabel="q (nm⁻¹)", ylabel="Intensity (a.u.)", ylabelfontsize=12, xlabelfontsize=12,
#         xtickfontsize=10, legendfontsize=10,
#           yticks=false, linewidth=2, leftmargin=3Plots.mm, xlim=(10, 50), grid=false)
# for i in axes(W, 2)
#     plot!(q, W[:,i], label="Basis $(i)", linewidth=2)
# end
# display(plt)

# ##### Plotting activations #####
# plt = plot(xlabel="Relative Location (μm)", ylabel="Activation (a.u.)", ylabelfontsize=12, xlabelfontsize=12,
#         xtickfontsize=10, ytickfontsize=10, legendfontsize=10,
#         linewidth=2, leftmargin=3Plots.mm, legend=:topright, grid=false)
# for i in axes(H, 1)
#     plot!(LinRange(-.75, .75, 76), H[i,77:152], label="Basis $(i)", linewidth=2)
# end
# display(plt)

phase_path = "paper/data/TaSnO/TaSnO_sticks.csv"
open(phase_path, "r") do f
    global cs = CrystalPhase(f, 0.1, FixedPseudoVoigt(0.1))
end

# optimize each W
# check amorphous
result_nodes = Vector{Node}(undef, RANK)

for i in 1:RANK
    # check amorphous
    bg = BackgroundModel(q, EQ(), 5, 100)
    W[:, i] ./= maximum(W[:,i])
    w = Wildcard([20., 38.], [.6, 0.2],  [2., 3.], "Amorphous", Lorentz(), [2., 2., 1., 1., .2, .5])
    amorphous = PhaseModel(nothing, nothing, bg)
    a = optimize!(amorphous, q, W[:, i], std_noise, mean_θ, std_θ, method=LM, objective="LS",
            maxiter=512, optimize_mode=Simple, regularization=true, verbose=false)


    res = norm(W[:, i] - evaluate!(zero(q), a, q))

    if res > 1.
        b = mcbl(W[:, i], q, 7)
        W[:, i] -= b
        W[W[:,i] .< 1e-5, i] .= 1e-5
        LT = Lazytree(cs, q)

        results = search!(LT, q, W[:, i], 2, 2, 1., false, false, .5, std_noise, mean_θ, std_θ,
                        method=LM, objective="LS", maxiter=128,
                        regularization=true)

        results = results[2:end]
        results = reduce(vcat, results)
        probs = get_probabilities(results, q, W[:, i], mean_θ, std_θ)
        result_node = results[argmax(probs)]

        result_nodes[i] = result_node

    end
end

phases = Vector{Set{String}}(undef, 201)

# Assignment of Phase related to space
H_THRESH = 0.2
FRAC_THRESH = 0.1
for colindex in 1:size(H, 2)
    phase_temp = Set{String}()
    for i in 1:size(H, 1)
        if H[i, colindex] > H_THRESH && isassigned(result_nodes, i)
            fractions = get_fraction(result_nodes[i].phase_model.CPs)
            for phase_idx in eachindex(fractions)
                if fractions[phase_idx] > FRAC_THRESH
                    push!(phase_temp, getproperty(result_nodes[i].phase_model.CPs[phase_idx], :name))
                end
            end
        end
    end
    phases[colindex] = phase_temp
end



# # Optimize each position with the identified phase only
pm_at_each_position = Vector{PhaseModel}(undef, 201)
# TODO: Parallelize?
bg = BackgroundModel(q, EQ(), 5, 100)
for i in tqdm(1:201)
    pm = get_phase_model_with_phase_names(phases[i], cs, bg)
    if !isempty(pm.CPs)
        normalized_data = data[336, i, 1:900] / maximum(data[336, i, 1:900])
        opt_pm = optimize!(pm, q, normalized_data,
                            std_noise, mean_θ, std_θ,
                            method=LM, objective="LS",
                            maxiter=512, regularization=true,
                            optimize_mode=Simple,
                            verbose=false)
        pm_at_each_position[i] = opt_pm
    end
end


# Plots
# All pattern with optimized phase
# for i in 1:201
#     plt = plot(q, data[336, i, 1:900]/ maximum(data[336, i, 1:900]), title="$(i)")
#     if isassigned(pm_at_each_position, i)
#         plot!(q, evaluate!(zero(q), pm_at_each_position[i], q))
#     end
#     display(plt)
# end

# Collect phases at their corresponding locations
# And plot their free lattice parameters
unique_phases = Set{String}()

plt = plot(xlabel="Relative Location (mm)", ylabel="Strain (%)", legend=:outerright, ylabelfontsize=12,
 xlabelfontsize=12, xtickfontsize=10, ytickfontsize=10, legendfontsize=10)

for i in unique(phases)
    if !isempty(i)
        push!(unique_phases, i...)
    end
end
simple_phase_name = Dict([("SnO2_P42/mnm", "SnO₂"),
                     ("Ta2O5_Pccm", "Ta₂O₅")])
for phase in unique_phases
    crystals = Vector{CrystalPhase}(undef, 201)
    for i in 1:201
        if isassigned(pm_at_each_position, i)
            exist_phases = getproperty.(pm_at_each_position[i].CPs, :name)
            ind = findall(x->x==phase, exist_phases)
            if !isempty(ind)
                crystals[i] = pm_at_each_position[i].CPs[ind[1]]
            end
        end
    end


    x = collect(1:201)[findall(isassigned.((crystals,), collect(1:201)))]
    println(x)
    free_param_num = crystals[x[1]].cl.free_param
    println(free_param_num)
    for i in 1:free_param_num
        param = [get_free_params(crystals[j])[i] for j in x]
        param_label = param_name[free_param_num][i]
        origin_param = getproperty(crystals[x[1]].origin_cl, eval(Meta.parse(":" * param_label)))
        plot!(0.01(x.-115), 100((param./origin_param).-1), label="$(simple_phase_name[phase]) $(param_label)", linewidth=2)
        # plot!(x, 100((param./maximum(param))), label="$(simple_phase_name[phase]) $(param_label)")
    end
end

display(plt)