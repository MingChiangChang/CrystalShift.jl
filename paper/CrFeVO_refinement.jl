using CSV
using DataFrames
using StatsBase

using CrystalShift
using CrystalShift: Lorentz, optimize!, PhaseModel, full_optimize!, FixedPseudoVoigt, PseudoVoigt
using CrystalShift: Wildcard, Crystal, twoθ2q, volume, get_lm_objective_func

using CovarianceFunctions: EQ, Matern
using Plots
using LaTeXStrings
using Measurements
using ForwardDiff
using LinearAlgebra
using LazyInverses
using ProgressBars

default(grid=false, fontfamily="Helvetica")

function make_string_to_numbers(ls)
    ls = ls[findall(x->x=='=', ls)[1]+1:end]
    ls = split(ls,',')
    return map(x->parse(Float64, x), map(String, ls))
end


d = ""
Q_IND = 22 # Q is in 22 row for the udi file
q_min = 200
q_max = 800
FeCr_ratio = [5.06757, 3.62887, 2.53125, 1.75, 1.42932, 1.16432, 1.02212, 0.79845, 0.673913, 0.422713, 0.28169]
sn = ["7134", "8743", "10205", "11853", "11845", "11836", "11828", "11819", "11811", "10121", "8616"]


# path = "/Users/ming/Downloads/CrFeV_toCornell/V54_xrd_A2theta_Int.csv"
# df = CSV.read(path, DataFrame)

phase_path = "paper/data/CrFeVO/sticks.csv"
open(phase_path, "r") do f
    global cs = CrystalPhase(f, 0.1, FixedPseudoVoigt(0.5))
end


#W = Wildcard([16.], [0.2],  [2.], "Amorphous", Lorentz(), [2.,  1., .5])


# comp = FeCr_ratio ./ (FeCr_ratio .+ 1)
Fe = (100 ./ (1 .+ FeCr_ratio) .* FeCr_ratio)
Fe_ = [i ± 5. for i in Fe]
Cr = 100 .- Fe
Cr_ = [i ± 5. for i in Cr]
comp = Fe_ ./ (Fe_ .+ Cr_)



##### Read and prepare data #####
# path = "paper/data/CrFeVO/ana__9_3594.udi"
# open(path) do f
#     global d = read(f, String)
# end
# split_data = split(d, '\n')
# sample_number = split(split_data[13][11:end-1], ',')
# q = make_string_to_numbers(split_data[Q_IND])[q_min:q_max]

##### Define containers #####
uncertainties = Vector{Vector{Float64}}()
cp = CrystalPhase[]
v = Float64[]

## CrystalShift settings ##
std_noise = .01
mean_θ = [1., .5, .1, 1., .5, .1]
std_θ = [.1, 10., .2, 0.05, 10., .2]
mod_peak_num = 200
loop_num = 32
peak_shift_iter = 256
peak_mod_iter=16
method = LM
objective = "LS"
regularization=true
peak_mod_mean=[1., 1.]
peak_mod_std=[.5, .05]

FeCr_ratio = [5.06757, 3.62887, 2.53125, 1.75, 1.42932, 1.16432, 1.02212, 0.79845, 0.673913, 0.422713, 0.28169]
q = npzread("paper/data/CrFeVO/CrFeVO_q.npy")
data = npzread("paper/data/CrFeVO/CrFeVO_data.npy")

for I in eachrow(data)

    # ind = findall(x->x==i, sample_number)[1]
    # I = make_string_to_numbers(split_data[Q_IND+ind])[q_min:q_max]
    # I ./= maximum(I)

    bg = BackgroundModel(q, Matern(3/2), 8., 10., rank_tol=1e-3)
    pm = PhaseModel(cs[1:2], nothing, bg)
    opt_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ,
                                            32, true,
                                            method, objective, Simple, 8, 1., false, 1e-6)

    @time opt_pm = full_optimize!(pm, q, I, std_noise, mean_θ, std_θ,
                            mod_peak_num = mod_peak_num,
                            loop_num=loop_num,
                            method=method,
                            peak_shift_iter=peak_shift_iter,
                            peak_mod_iter=peak_mod_iter,
                            objective=objective,
                            regularization=regularization,
                            peak_mod_mean=peak_mod_mean,
                            peak_mod_std=peak_mod_std)


    # Plot fitting spectrum
    plt = plot(q, I, label="Raw Data", legend=:topleft, xlim=(16.5, 32.535), ylim=(0.0, 1.05), linewidth=10,
                xlabel="q (nm⁻¹)", ylabel="Normalized Intensitiy", labelfontsize=20, xtickfontsize=16, ytickfontsize=16, titlefontsize=20, legendfontsize=16)
    # savefig("raw_data.png")
    plot!(q, evaluate!(zero(q), opt_pm.background, q), label="Background", linewidth=4)
    # savefig("raw_bg.png")
    plot!(q, evaluate!(zero(q), opt_pm.CPs[2], q), label="SnO₂", linewidth=6)
    # savefig("with_bg_sno.png")
    plot!(q, evaluate!(zero(q), opt_pm.CPs[1], q), label="CrₓFe₍₀.₅₋ₓ₎VO₄", linewidth=6)
    plot!(q, evaluate!(zero(q), opt_pm, q), color=:red, label="Optimized Result", linewidth=4)

    plot!(size=(900,900), left_margin=5Plots.mm, bottom_margin=5Plots.mm, dpi=300, framestyle = :box, legend=:topright)
    # savefig("final.svg")
    display(plt)

    ##### Uncertainty evaluateion #####
    I = I - evaluate!(zero(q), opt_pm.background, q)
    opt_pm = PhaseModel(opt_pm.CPs, nothing, nothing)
    params = get_free_params(pm)
    f = get_lm_objective_func(opt_pm, q, I, zero(I), opt_stn)
    r = zeros(Real, length(I) + length(params))
    function res(log_θ)
        sum(abs2, f(r, log_θ)/maximum(I))
    end
    log_θ = log.(get_free_params(opt_pm))

    H = ForwardDiff.hessian(res, log_θ)
    val = res(log_θ)
    uncer = sqrt.(diag(val / (length(q) - length(log_θ)) * inverse(H)))

    push!(uncertainties, uncer)
    push!(cp, opt_pm.CPs[1])
    push!(v, volume(opt_pm.CPs[1].cl))
end

c_a = [exp(log(cp[i].cl.a) ± 2uncertainties[i][1]) for i in eachindex(cp)]
norm_a = cs[1].cl.a
norm_c_a = (c_a .-norm_a) ./norm_a .* 100
c_b = [exp(log(cp[i].cl.b) ± 2uncertainties[i][2]) for i in eachindex(cp)]
norm_b = cs[1].cl.b
norm_c_b = (c_b .-norm_b) ./norm_b .* 100
c_c = [exp(log(cp[i].cl.c) ± 2uncertainties[i][3]) for i in eachindex(cp)]
norm_c = cs[1].cl.c
norm_c_c = (c_c .-norm_c) ./norm_c .* 100
c_β = [exp(log(cp[i].cl.β) ± 2uncertainties[i][4]) for i in eachindex(cp)]
norm_β = cs[1].cl.β
norm_c_β = (c_β.-norm_β) ./norm_β .* 100
cl = [norm_c_a, norm_c_b, norm_c_c, norm_c_β]
lattice_parm = ["a", "b", "c", "β"]
plt = plot(layout = 4, legend=false)

xticks = [1,2,3,4,5,6]
default(labelfontsize=20, xtickfontsize=16, ytickfontsize=16, titlefontsize=20,linewidth=6)
x_val = getproperty.(comp, :val)
x_err = getproperty.(comp, :err)
p1 = plot(comp, cl[1], xerr=x_err, yerr=getproperty.(c_a, :err)./norm_a.*100, title="a", xlabel="", ylabel="Lattice Strain (%)", color=:red,ylims=(-1.0, 1.0))
scatter!(comp, cl[1], color=:pink)
p2 = plot(comp, cl[2], xerr=x_err, yerr=getproperty.(c_b, :err)./norm_b.*100, title="b", xlabel="", ylabel="Lattice Strain (%)", color=:orange, ylims=(-1.0, 1.0))
scatter!(comp, cl[2], color=:gold)
p3 = plot(comp, cl[3], xerr=x_err, yerr=getproperty.(c_c, :err)./norm_c.*100, title="c", ylabel="Lattice Strain (%)", xlabel="Fe/(Fe+Cr)", color=:green, ylims=(-1.0, 1.0))
scatter!(comp, cl[3], color=:lightgreen)
p4 = plot(comp, getproperty.(c_β, :val).* 180 ./pi, xerr=x_err, yerr=getproperty.(c_β, :err).* 180 ./pi, title="β", xlabel="Fe/(Fe+Cr)",ylabel="β (°)", color=:blue,
       ylim=(mean(getproperty.(c_β, :val).* 180 ./pi)*0.99, mean(getproperty.(c_β, :val).* 180 ./pi)*1.01 ), ytickfontsize=12)
scatter!(comp, getproperty.(c_β, :val).* 180 ./pi, color=:cornflowerblue)
plt = plot(p1, p2, p3, p4, layout = (2, 2), legend=false, left_margin=5Plots.mm, bottom_margin=2Plots.mm)
plot!(size=(800,1000), dpi=300, framestyle = :box ,xticks=[0.2, 0.4, 0.6, 0.8], xlims=(0.2, 0.85))
# savefig("CrFeVO_lattice")
display(plt)

uncer_v = Vector{Measurement}()
for i in 1:11
    push!(uncer_v, volume(Monoclinic{Measurement}(c_a[i], c_b[i], c_c[i], c_β[i])))
end

xtick=collect(1:5)
plot(comp, uncer_v, xlabel="Fe/(Fe+Cr)", ylabel="Unit Cell volume (Å³)", legend=false)
scatter!(comp, v, color=:blue)



## CIF sticks
function rectangle(pos, w, h)
    x, y = pos
    Plots.Shape([x-w/2, x+w/2, x+w/2, x-w/2], [y, y, y+h, y+h])
end

phase_path = "paper/data/CrFeVO/sticks.csv"
open(phase_path, "r") do f
    global cs = CrystalPhase(f, 0.1, FixedPseudoVoigt(0.5))
end

plt = plot(size=(900, 200), xlim=(16.5, 32.535), ylim=(0, 1), xtickfontsize=16,
          legendfontsize=12, frame=:box, bottom_margin=10Plots.mm, left_margin=26.8Plots.mm, yticks=false, labelfontsize=20, ylabel="(a.u)", xlabel="q (nm⁻¹)")
names = ["SnO₂ 01-070-6153", "Cr₀.₅Fe₀.₅VO₄ 04-011-4573"]
rcs = reverse(cs)
for i in eachindex(rcs)
    peak_qs = (rcs[i].cl).(rcs[i].peaks).*10
    # plt = plot()
    for j in eachindex(peak_qs)
        if j == 1
            plot!(rectangle([peak_qs[j], 0.01], 0.1, rcs[i].peaks[j].I), c=palette(:auto)[i+2], linewidth=0, label=names[i])
        end
        plot!(rectangle([peak_qs[j], 0.001], 0.1, rcs[i].peaks[j].I), c=palette(:auto)[i+2], linewidth=0, label=false)
    end
end
# savefig("reference_sticks.svg")
display(plt)

lattice_a = CSV.read("paper/data/CrFeVO/reference_lp_a.csv", DataFrame)
lattice_b = CSV.read("paper/data/CrFeVO/reference_lp_b.csv", DataFrame)
lattice_c = CSV.read("paper/data/CrFeVO/reference_lp_c.csv", DataFrame)
a = lattice_a[!, 2] ./= 9.8245
b = lattice_b[!, 2] ./= 8.8776
c = lattice_c[!, 2] ./= 6.8252

a .-= 1
b .-= 1
c .-= 1
a .*= 100
b .*= 100
c .*= 100
a = reverse(a)
b = reverse(b)
c = reverse(c)

default(labelfontsize=20, xtickfontsize=16, ytickfontsize=16, titlefontsize=20,linewidth=6)
x_val = getproperty.(comp, :val)
x_err = getproperty.(comp, :err)
p1 = plot(comp, cl[1], xerr=x_err, yerr=getproperty.(c_a, :err)./norm_a.*100, title="a", xlabel="Fe/(Fe+Cr)", ylabel="Lattice Strain (%)", color=:red,ylims=(-1.0, 1.0), label="CrystalShift")
scatter!(comp, cl[1], color=:pink, label=nothing)
plot!(comp, a, color=:pink, alpha=0.5, label="Zhou et al.")
scatter!(comp, a, color=:pink, alpha=0.5, label=nothing)
p2 = plot(comp, cl[2], xerr=x_err, yerr=getproperty.(c_b, :err)./norm_b.*100, title="b", xlabel="Fe/(Fe+Cr)", ylabel="Lattice Strain (%)", color=:orange, ylims=(-1.0, 1.0), label="CrystalShift")
scatter!(comp, cl[2], color=:gold, label=nothing)
plot!(comp, b, color=:gold, markershape=:circle, alpha=0.5, label="Zhou et al.")
scatter!(comp, b, color=:gold, label=nothing)
p3 = plot(comp, cl[3], xerr=x_err, yerr=getproperty.(c_c, :err)./norm_c.*100, title="c", ylabel="Lattice Strain (%)", xlabel="Fe/(Fe+Cr)", color=:green, ylims=(-1.0, 1.0), label="CrystalShift")
scatter!(comp, cl[3], color=:lightgreen, label=nothing)
plot!(comp, c, color=:lightgreen, markershape=:circle, alpha=0.5, label="Zhou et al.")
scatter!(comp, c, color=:lightgreen, label=nothing)
p4 = plot(comp, getproperty.(c_β, :val).* 180 ./pi, xerr=x_err, yerr=getproperty.(c_β, :err).* 180 ./pi, title="β", xlabel="Fe/(Fe+Cr)",ylabel="β (°)", color=:blue,
       ylim=(mean(getproperty.(c_β, :val).* 180 ./pi)*0.99, mean(getproperty.(c_β, :val).* 180 ./pi)*1.01 ), ytickfontsize=12)
scatter!(comp, getproperty.(c_β, :val).* 180 ./pi, color=:cornflowerblue)
plt = plot(p1, p2, p3, layout = (1, 3), legend=false, left_margin=5Plots.mm, bottom_margin=10Plots.mm)
plot!(size=(1200,600), dpi=300, framestyle = :box ,xticks=[0.2, 0.4, 0.6, 0.8], xlims=(0.2, 0.85), legend=true, legendfontsize=12)
# savefig("CrFeVO_lattice")
display(plt)