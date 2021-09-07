using Plots
using JSON
using ProgressBars

include("../src/CrystalShift.jl")

function get_unique_phase_names(wafer_result)
    names = String[]
    println(typeof(wafer_result))
    for stripe in wafer_result
        for phase in stripe["phase_results"]
            n = phase["phase_name"]
            if !(n in names) & phase["isCenter"]
                push!(names, n)
            end
        end
    end
    names
end

function get_center(stripe_dict)
    for phase in stripe_dict["phase_results"]
        if phase["isCenter"]
            return phase
        end
    end
end

path = "/Users/r2121/Desktop/Code/Crystallography_based_shifting/data/"
path = path * "TaSnO.json"

results = JSON.parsefile(path)
names = get_unique_phase_names(results)

plt = plot()

for name in names
    tpeak = Float64[]
    dwell = Float64[]
    cation = Float64[]
    for r in tqdm(results)
        center_result = get_center(r)
        if !isnothing(center_result) && center_result["phase_name"] == name
            println(typeof(r["tpeak"]), typeof(r["dwell"]), typeof(r["cation_ratio"]), r["cation_ratio"])
            push!(tpeak, r["tpeak"])
            push!(dwell, r["dwell"])
            push!(cation, r["cation_ratio"][1])
        end
    end
    scatter!(cation, tpeak, label=name)
end

display(plt)
