using CrystalShift
# using DelimitedFiles
using LinearAlgebra

using CrystalShift: CrystalPhase, optimize!, Lorentz, PseudoVoigt, FixedPseudoVoigt, get_free_lattice_params
using CrystalShift: evaluate

using Test
using Plots

function synthesize_multiphase_data(cps::AbstractVector{<:CrystalPhase},
                                    x::AbstractVector)
    r = zero(x)
    full_params = Float64[]
    interval_size = 0.01

    for cp in cps
        params = get_free_lattice_params(cp)

        scaling = (interval_size.*rand(size(params, 1)).-interval_size/2).+1
        @. params = params*scaling
        params = vcat(params, 0.5.+3rand(1), 0.1.+0.1(rand(1)))
        full_params = vcat(full_params, params)
    end
    r = evaluate(cps, full_params, x)

    r/maximum(r), full_params
end

function test_multiphase_optimize(cps::AbstractVector{<:CrystalPhase},
                                x::AbstractVector, num_phase::Int,
                                optimize_mode::OptimizationMode,
                                method::OptimizationMethods, maxiter::Int=64, objective::String = "LS",
                                verbose = false)
    phase = rand(1:size(cps, 1), num_phase)
    y, full_params = synthesize_multiphase_data(cps[phase], x)
    noise = 0.03rand(size(x, 1))
    y += noise
    y /= maximum(y)
    c = optimize!(PhaseModel(cps[phase]), x, y, std_noise, mean_θ, std_θ;
                objective = objective, optimize_mode = optimize_mode,
                method = method, maxiter = maxiter, em_loop_num=2,
                regularization = true, verbose = verbose)
    #for i in c
    #    println(get_free_params(i))
    #end
    #println(full_params)
    if isa(c, Tuple)
        c, _ = c
    end
    plt = plot(x, y)
    plot!(x, evaluate!(zero(x), c, x))
    display(plt)
    norm((c.CPs).(x).-y)
end

std_noise = .05
mean_θ = [1., 1., .2]
std_θ = [.5, .5, 1.]

# CrystalPhas object creation
path = "data/"
phase_path = path * "sticks.csv"
phase_path = "/Users/ming/Downloads/YourCustomFileName/sticks.csv"
f = open(phase_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n") # Windows: #\r\n ...
else
    s = split(read(f, String), "#\n")
end

if s[end] == ""
    pop!(s)
end

cs = Vector{CrystalPhase}(undef, size(s))
@. cs = CrystalPhase(String(s), (0.1, ), (FixedPseudoVoigt(0.5), )) # For ease of testing fast

x = collect(10.:.1:80.)
y = zero(x)
evaluate!(y, cs[1], x)
evaluate!(y, cs[2], x)

test_multiphase_optimize(cs[[3,5]], x, 3, EM, LM)


# LT = Lazytree(cs, x, 3)#, 20,true)

# @time phase = search!(LT, x, y, 2, 3, std_noise, mean_θ, std_θ,
#                     maxiter=128, regularization=true)
# println("")