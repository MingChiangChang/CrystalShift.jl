using Plots
using JSON
using ProgressBars

using CrystalShift
using CrystalShift: Lorentz

function get_unique_phase_names(wafer_result)
    names = String[]

    for stripe in wafer_result
        for basis in stripe["phase_results"]
            for phase in basis
                n = phase["phase_name"]
                if !(n in names) & phase["isCenter"]
                    push!(names, n)
                end
            end
        end
    end
    names
end

function get_center(stripe_dict)
    for basis in stripe_dict["phase_results"]
        for phase in basis
            if phase["isCenter"]
                return basis
            end
        end
    end
end

function get_phase_names(phase_results)
    [pr["phase_name"] for pr in phase_results]
end

function get_lattice_params(phase_result, name)
    for pr in phase_result
        if pr["phase_name"] == name
            phase = pr["phase"][1]["cl"]
            return [phase["a"], phase["b"], phase["c"], phase["α"], phase["β"], phase["γ"]]
        end
    end
end

# Fill in json path from running TaSnO.jl
path = "data/TaSnO_metric_2023-11-23_11:17.json"
results = JSON.parsefile(path)
names = get_unique_phase_names(results)

plt = plot(legend=(0.65, 0.04), xlims=(0, 1), fontsize=10, size=(600, 400))
phase_names = ["SnO2_P42/mnm", "Ta2O5_Pccm", "Pyrochlore_Fd-3mZ"]
chem_names = ["SnO₂_P4₂/mnm", "Ta₂O₅_Pccm", "Ta₂Sn₂O₇_Fd-3mZ"] # For legends

s = [12, 9, 7]
color = [10, 8, 6]

for (idx, phase_name) in enumerate(phase_names)
    global tpeak = Float64[]
    global dwell = Float64[]
    global cation = Float64[]
    global lattice = Vector{Vector{Float64}}()
    for (idx, r) in enumerate(results)
        center_result = get_center(r)

        if !isnothing(center_result)
            result_names = get_phase_names(center_result)
            if phase_name in result_names
                #println("$(idx) $(r["tpeak"]) $(r["dwell"])")
                push!(lattice, get_lattice_params(center_result, phase_name))
                push!(tpeak, r["tpeak"])
                push!(dwell, r["dwell"])
                push!(cation, r["cation_ratio"][1])
            end
        end
    end
    scatter!(cation, tpeak, label=chem_names[idx], markersize=color[idx], legend=true, frame=:box,
             xtickfontsize=12, xlabelfontsize=11, ylabelfontsize=11, ytickfontsize=10, right_margin=3Plots.mm)

    ylabel!("T peak (°C)")
    xlabel!("Ta/(Ta+Sn)")
end
display(plt)
